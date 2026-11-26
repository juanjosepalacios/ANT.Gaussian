// Reveal a keyword when following a link into a collapsed reference entry.
function revealLinkedEntry() {
  const identifier = decodeURIComponent(window.location.hash.slice(1));
  const target = document.getElementById(identifier);
  if (!target) return;
  let entry = target.closest('details');
  while (entry) {
    entry.open = true;
    entry = entry.parentElement.closest('details');
  }
  target.scrollIntoView({ block: 'start' });
}
window.addEventListener('hashchange', revealLinkedEntry);
window.addEventListener('DOMContentLoaded', revealLinkedEntry);

window.addEventListener('DOMContentLoaded', () => {
  const header = document.querySelector('.section-header');
  const updateHeaderHeight = () => {
    document.documentElement.style.setProperty('--section-nav-height', `${header.offsetHeight}px`);
  };
  updateHeaderHeight();
  new ResizeObserver(updateHeaderHeight).observe(header);
});

// Keep the selected Gaussian interface when navigating shared documentation.
window.addEventListener('DOMContentLoaded', () => {
  if (new URLSearchParams(window.location.search).get('code') !== 'gaussian16') return;
  const destinations = {
    'ant_g03.html': 'ant_g03.html?code=gaussian16',
    'manual_1_g03.html': 'g16-installation.html?v=2',
    'manual_3_g03.html': 'manual_3_g03.html?code=gaussian16',
    'input_files_g03.html': 'g16-input-files.html?v=2',
    'publications_g03.html': 'publications_g03.html?code=gaussian16',
    'examples_g03.html': 'g16-examples.html'
  };
  document.querySelectorAll('.section-header nav a, .pager a').forEach(link => {
    const href = link.getAttribute('href');
    if (destinations[href]) link.setAttribute('href', destinations[href]);
  });
  const nav = document.querySelector('.section-header nav');
  if (nav) nav.setAttribute('aria-label', 'ANT.Gaussian16 sections');
});
