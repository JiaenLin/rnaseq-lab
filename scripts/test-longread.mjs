// Regression tests for the isoform layer, against the exact header
// wf-transcriptomes writes to out/cohort/transcript_counts.tsv.
import { parseMatrix } from '../src/lib/matrix.ts'
import {
  isTranscriptMatrix, readLongRead, applySqanti, displayName, transcriptsCsv, describe,
} from '../src/lib/longread.ts'

let failed = 0
const check = (name, got, want) => {
  const ok = JSON.stringify(got) === JSON.stringify(want)
  if (!ok) failed++
  console.log(`${ok ? '  ok  ' : '  FAIL'} ${name}${ok ? '' : `\n        got  ${JSON.stringify(got)}\n        want ${JSON.stringify(want)}`}`)
}
const rows = (...lines) => lines.map(l => l.split('\t'))

const HEAD = 'TXNAME\tGENEID\tNDR\tnovelGene\tnovelTranscript\ttxClassDescription\treadCount\trelReadCount\trelSubsetCount\ttxid\teqClassById\tgene_name\ttranscript_name\tA24M_4\tA24M_5\tWT_1\tWT_2'
const TAB = rows(
  HEAD,
  // annotated, two isoforms of one gene
  'ENSMUST00000000001\tENSMUSG00000000001\tNA\tFALSE\tFALSE\tannotation\t900\t1\t1\t1\t1\tGnai3\tGnai3-201\t10\t12\t40\t44',
  'ENSMUST00000000002\tENSMUSG00000000001\tNA\tFALSE\tFALSE\tannotation\t400\t1\t1\t2\t2\tGnai3\tGnai3-202\t30\t28\t5\t6',
  // novel, same gene
  'BambuTx7\tENSMUSG00000000001\t0.01\tFALSE\tTRUE\tnewWithin\t41\t0.74\t1\t3\t3\tGnai3\tNA\t3\t4\t1\t0',
  // novel gene: no symbol at all
  'BambuTx9\tBambuGene103\t0.02\tTRUE\tTRUE\tnewGene-spliced\t9\t0.5\t1\t4\t4\tNA\tNA\t2\t1\t0\t1',
)

console.log('\nRECOGNISING A TRANSCRIPT MATRIX')
check('TXNAME is transcript-level', isTranscriptMatrix(TAB[0]), true)
check('gene_id is not', isTranscriptMatrix(['gene_id', 'gene_name', 'a', 'b']), false)
check('GENEID is not', isTranscriptMatrix(['GENEID', 'newGeneClass', 'a', 'b']), false)

const m = parseMatrix(TAB)
const lr = readLongRead(m, TAB)

console.log('\nREADABLE NAMES')
// The collision that shipped once: `tx_name` normalises to `txname`, which is
// TXNAME — the ID column — so every name became its own accession.
check('annotated transcript uses its reference name',
  lr.transcripts[0].display_name, 'Gnai3-201')
check('and is not its accession',
  lr.transcripts.filter(t => t.display_name === t.transcript_id).length, 1)
check('novel isoform of a known gene is named from the symbol',
  lr.transcripts[2].display_name, 'Gnai3-novel-1')
check('novel gene falls back to the accession',
  lr.transcripts[3].display_name, 'BambuTx9')
check('novel counter is per gene, not global',
  displayName('BambuTx99', 'Ttn', '', 1), 'Ttn-novel-1')

console.log('\nGENE MATRIX SUMMED FROM THE TRANSCRIPT MATRIX')
check('genes, not transcripts', lr.nGenes, 2)
const gl = lr.geneCountsCsv.trim().split('\n')
check('header is gene_id + samples', gl[0], 'gene_id,A24M_4,A24M_5,WT_1,WT_2')
// 10+30+3 / 12+28+4 / 40+5+1 / 44+6+0
check('three isoforms summed into their gene', gl[1], 'ENSMUSG00000000001,43,44,46,50')
check('the novel gene stands alone', gl[2], 'BambuGene103,2,1,0,1')
check('symbols carried to the gene layer', lr.geneNames.get('ENSMUSG00000000001'), 'Gnai3')
check('every transcript maps to a gene', lr.txToGene.size, 4)

console.log('\nSQANTI JOIN')
const SQ = rows(
  'isoform\tchrom\tstructural_category\tassociated_gene',
  'ENSMUST00000000001\t1\tfull-splice_match\tGnai3',
  'BambuTx7\t1\tnovel_in_catalog\tGnai3',
)
const j = applySqanti(lr.transcripts, SQ)
check('only the rows SQANTI names are touched', j.matched, 2)
check('category replaced', lr.transcripts[0].structural_category, 'full-splice_match')
check('novel one too', lr.transcripts[2].structural_category, 'novel_in_catalog')
check('unmentioned row keeps what bambu said',
  lr.transcripts[1].structural_category, 'annotation')

let threw = ''
try { applySqanti(lr.transcripts, rows('a\tb', '1\t2')) } catch (e) { threw = e.message }
check('a file that is not a classification says so', /SQANTI3 classification/.test(threw), true)

console.log('\nBUNDLE FILE')
const csv = transcriptsCsv(lr.transcripts).trim().split('\n')
check('header', csv[0],
  'transcript_id,gene_id,gene_name,transcript_name,display_name,structural_category,novel')
check('one row per transcript, nothing merged', csv.length, 5)
check('novel flag is written', csv[3].endsWith(',TRUE'), true)

let bad = ''
try {
  readLongRead(parseMatrix(rows('TXNAME\ta\tb', 'T1\t1\t2', 'T2\t3\t4')),
    rows('TXNAME\ta\tb', 'T1\t1\t2', 'T2\t3\t4'))
} catch (e) { bad = e.message }
check('a transcript matrix with no gene column is refused', /needs a gene column/.test(bad), true)

console.log('\nDEFECTS FOUND IN REVIEW, PINNED')
{
  // A gene id carrying a comma is quoted by matrix.ts's writer. A naive
  // split(',') then shifted every sample one column left, on that row only.
  const T = rows(
    'TXNAME\tGENEID\tgene_name\ttranscript_name\ta\tb',
    'T1\tG,1\tGx\tGx-201\t10\t12',
    'T2\tG,1\tGx\tGx-202\t30\t28',
    'T3\tG2\tGy\tGy-201\t7\t7',
  )
  const l = readLongRead(parseMatrix(T), T)
  const g = l.geneCountsCsv.trim().split('\n')
  check('a comma in a gene id does not shift the counts', g[1], '"G,1",40,40')
  check('and the next gene is untouched', g[2], 'G2,7,7')
}
{
  // The transcript matrix shipped in the bundle must be keyed by transcript_id
  // and rounded — not matrix.ts's `gene_id`-headed, still-fractional CSV.
  const T = rows(
    'TXNAME\tGENEID\tgene_name\ttranscript_name\ta\tb',
    'T1\tG1\tGx\tGx-201\t36.69048\t20.83751',
    'T2\tG1\tGx\tGx-202\t0.4\t0.6',
  )
  const l = readLongRead(parseMatrix(T), T)
  const tx = l.txCountsCsv.trim().split('\n')
  check('transcript matrix is keyed by transcript_id', tx[0], 'transcript_id,a,b')
  check('and is rounded', tx[1], 'T1,37,21')
  // 37+0 and 21+1: the gene total is the sum of the SAME integers the
  // transcript file ships and the usage test is handed.
  check('the gene total is the sum of exactly those integers',
    l.geneCountsCsv.trim().split('\n')[1], 'G1,37,22')
}
{
  // structural_category speaks two incompatible vocabularies; the bundle has
  // to say which, and a partial SQANTI join is honestly 'mixed'.
  const T = rows(
    'TXNAME\tGENEID\ttxClassDescription\tgene_name\ttranscript_name\ta\tb',
    'T1\tG1\tnewWithin\tGx\tNA\t5\t6',
    'T2\tG1\tannotation\tGx\tGx-201\t5\t6',
  )
  const l = readLongRead(parseMatrix(T), T)
  check('bambu class strings are labelled bambu, not SQANTI', l.vocabulary, 'bambu')
  check('and describe() says so', /bambu classes/.test(describe(l)), true)
  const v = applySqanti(l.transcripts, rows('isoform\tstructural_category', 'T1\tnovel_in_catalog'))
  check('a join that reached only some rows is mixed', v.vocabulary, 'mixed')
  const l2 = readLongRead(parseMatrix(T), T)
  const v2 = applySqanti(l2.transcripts,
    rows('isoform\tstructural_category', 'T1\tnovel_in_catalog', 'T2\tfull-splice_match'))
  check('a complete join is sqanti', v2.vocabulary, 'sqanti')
}

console.log(failed ? `\n${failed} test(s) failed\n` : '\nAll isoform-layer tests passed\n')
process.exit(failed ? 1 : 0)
