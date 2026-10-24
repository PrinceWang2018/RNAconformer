// Desktop-shell bridge (pywebview): native Open/Save dialogs. Without the bridge the app falls back to the
// browser's file chooser, the File System Access save dialog, or downloads.
export class Native {
  constructor(api) {
    this.api = api;
    this.ready = !!(window.pywebview && window.pywebview.api);
    window.addEventListener('pywebviewready', () => { this.ready = true; });
  }

  available() { return this.ready && !!(window.pywebview && window.pywebview.api); }

  async save(blob, name, description = 'File') {
    const path = await window.pywebview.api.choose_save_path(name, description);
    if (!path) return null;
    const r = await fetch('/api/write-file', {
      method: 'POST', body: blob,
      headers: { 'X-Viewer-Token': this.api.token || '', 'X-Path': encodeURIComponent(path), 'Content-Type': 'application/octet-stream' },
    });
    if (!r.ok) throw new Error((await r.json().catch(() => ({}))).error || `could not write ${path}`);
    return path;
  }

  async chooseFiles(kind = 'structures') { return (await window.pywebview.api.choose_open_files(kind)) || []; }
  async readText(path) { return window.pywebview.api.read_text_file(path); }
  async chooseFolder() { return await window.pywebview.api.choose_folder(); }
  async reveal(path) { return window.pywebview.api.reveal(path); }
}
