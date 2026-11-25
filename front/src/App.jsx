import { useState } from 'react';
import Header from './componentes/Header';
import Home from './componentes/Home';
import About from './componentes/About';
import './App.css';

function App() {
  const [paginaAtual, setPaginaAtual] = useState('home');

  return (
    // Exactly one screen tall: header and footer stay put, and only <main> can
    // scroll (the About page). Home sizes itself to fit inside it.
    <div className="flex flex-col h-screen font-sans bg-white">
      {/* O Header recebe a função para mudar de página */}
      <div className="shrink-0">
        <Header paginaAtual={paginaAtual} setPaginaAtual={setPaginaAtual} />
      </div>
      
      <main className="flex-1 min-h-0 overflow-y-auto">
        {/* Renderização Condicional */}
        {paginaAtual === 'home' && <Home />}
        {paginaAtual === 'about' && <About />}
      </main>

      {/* Footer (Opcional, vindo do modelo do Stitch) */}
      <footer className="shrink-0 w-full bg-white border-t border-zinc-100 px-12 py-3">
        <div className="text-[10px] font-medium uppercase tracking-[0.2em] text-[#596061]">
          <a href="https://mcblab.github.io/" target="_blank" rel="noopener noreferrer" className="hover:text-[#2d3435] transition-colors">
            Made by MCB Lab at UFRN
          </a>
        </div>
      </footer>
    </div>
  );
}

export default App;