// Test harness for tools/viewer_smoketest.py: loads RCV in an iframe and runs one scenario from steps.js.
const q = new URLSearchParams(location.search);
const logEl = document.getElementById('log');
let chain = Promise.resolve();       // keep log lines in order (DONE must arrive last)
const log = (...a) => {
  logEl.textContent += a.join(' ') + '\n';
  const body = JSON.stringify({ level: 'test', message: `[${q.get('steps')}] ` + a.join(' ') });
  chain = chain.then(() => fetch('/api/client-log', { method: 'POST', headers: { 'X-Viewer-Token': q.get('token'), 'Content-Type': 'application/json' }, body }).catch(() => {}));
};
const f = document.getElementById('f');
f.src = '/' + (q.get('lang') ? '?lang=' + encodeURIComponent(q.get('lang')) : '') + '#token=' + q.get('token');
const sleep = ms => new Promise(r => setTimeout(r, ms));
(async () => {
  try {
    for (let i = 0; i < 200; i++) { if (f.contentWindow?.rcv && f.contentDocument?.body?.dataset.ready) break; await sleep(100); }
    const app = f.contentWindow.rcv;
    f.contentWindow.__noFilePicker = true;
    const steps = await import('/_dev/steps.js');
    await steps[q.get('steps')](app, f.contentWindow, log, sleep);
    log('DONE');
  } catch (e) { log('ERROR', e.message, e.stack); }
  document.body.dataset.done = '1';
})();
