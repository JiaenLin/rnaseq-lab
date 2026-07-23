import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

// base './' → works under GitHub Pages /<repo>/. COOP/COEP headers in dev give
// SharedArrayBuffer locally (on Pages the coi-serviceworker provides it).
export default defineConfig({
  base: './',
  plugins: [react()],
  server: {
    headers: {
      'Cross-Origin-Opener-Policy': 'same-origin',
      'Cross-Origin-Embedder-Policy': 'require-corp',
    },
  },
})
