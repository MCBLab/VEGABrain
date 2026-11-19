import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

// https://vite.dev/config/
export default defineConfig({
  plugins: [react()],
  // O site e publicado em https://mcblab.github.io/VEGABrain/, entao os assets
  // precisam ser referenciados a partir de /VEGABrain/ e nao da raiz. Se o
  // repositorio for renomeado de novo, este caminho tem que acompanhar.
  base: '/VEGABrain/',
})
