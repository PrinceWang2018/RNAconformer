"""Check for updates (rcv/update.py): versions, the installer for each platform, the website fallback, auto-check."""
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from rcv import update as U

ROOT = Path(__file__).resolve().parents[1]
NAMES = ['RCV-0.8.0-Linux-amd64.deb', 'RCV-0.8.0-Linux-x86_64.rpm', 'RCV-0.8.0-macOS-arm64.dmg',
         'RCV-0.8.0-macOS-x86_64.dmg', 'RCV-0.8.0-Windows-x64-Setup.exe', 'SHA256SUMS.txt']
GITHUB = {'tag_name': 'v0.8.0', 'published_at': '2026-12-01T10:00:00Z',
          'html_url': 'https://github.com/x/y/releases/tag/v0.8.0', 'assets': [{'name': n} for n in NAMES]}
STATS = {'downloads': {'releases': [
    {'repo': U.RELEASE_REPO, 'tag': 'v0.7.10', 'date': '2026-11-20', 'assets': [{'name': 'RCV-0.7.10-Linux-x86_64.rpm'}]},
    {'repo': U.RELEASE_REPO, 'tag': 'v0.7.9', 'date': '2026-11-21', 'assets': []},
    {'repo': 'PrinceWang2018/RNAconformer', 'tag': 'v9.9.9', 'date': '2026-11-22', 'assets': []}]}}


def fake(github=GITHUB, stats=STATS):
    def fetch(url):
        data = github if 'api.github.com' in url else stats
        if isinstance(data, Exception):
            raise data
        return data
    return fetch


class TestVersions(unittest.TestCase):
    def test_compare(self):
        self.assertEqual(U.parse_version('v0.7.10'), (0, 7, 10))
        self.assertTrue(U.is_newer('0.7.10', '0.7.9'))          # numbers, not text
        self.assertTrue(U.is_newer('v0.8', '0.7.4'))
        self.assertFalse(U.is_newer('0.7.3', '0.7.4'))          # a development version ahead of the release
        self.assertFalse(U.is_newer('0.7.4', '0.7.4.0'))
        self.assertFalse(U.is_newer('', '0.7.4'))


class TestPlatform(unittest.TestCase):
    def test_keys(self):
        self.assertEqual(U.platform_key('win32', 'AMD64'), 'windows-x64')
        self.assertEqual(U.platform_key('darwin', 'arm64'), 'macos-arm64')
        with mock.patch.object(U, '_rosetta', return_value=False):
            self.assertEqual(U.platform_key('darwin', 'x86_64'), 'macos-x86_64')
        with mock.patch.object(U, '_rosetta', return_value=True):     # Intel build on an Apple Silicon Mac
            self.assertEqual(U.platform_key('darwin', 'x86_64'), 'macos-arm64')
        self.assertEqual(U.platform_key('linux', 'x86_64', {'ID': 'ubuntu', 'ID_LIKE': 'debian'}, False), 'linux-deb')
        self.assertEqual(U.platform_key('linux', 'x86_64', {'ID': 'rocky', 'ID_LIKE': 'rhel centos fedora'}, False), 'linux-rpm')
        self.assertEqual(U.platform_key('linux', 'x86_64', {'ID': 'opensuse-leap', 'ID_LIKE': 'suse opensuse'}, False), 'linux-rpm')
        self.assertEqual(U.platform_key('linux', 'x86_64', {'ID': 'fedora'}, True), 'linux-deb')   # installed from the .deb
        self.assertIsNone(U.platform_key('linux', 'aarch64', {'ID': 'ubuntu'}, False))           # no ARM Linux build

    def test_every_platform_finds_its_installer(self):
        for key in U.PLATFORMS:
            r = U.check('0.7.4', fetch=fake(), key=key, installed=True)
            self.assertTrue(r['newer'])
            self.assertEqual(r['latest'], '0.8.0')
            self.assertIsNotNone(r['file'], key)
            self.assertEqual(r['download'], f"{U.WEBSITE}download/v0.8.0/{r['file']}")
        self.assertEqual(U.check('0.7.4', fetch=fake(), key='linux-rpm', installed=True)['file'],
                         'RCV-0.8.0-Linux-x86_64.rpm')
        none = U.check('0.7.4', fetch=fake(), key=None, installed=True)
        self.assertIsNone(none['download'])
        self.assertTrue(none['page'].endswith('#download'))


class TestSources(unittest.TestCase):
    def test_github_first(self):
        r = U.latest_release(fake())
        self.assertEqual((r['version'], r['source'], r['date']), ('0.8.0', 'github', '2026-12-01'))

    def test_website_when_github_fails(self):
        r = U.latest_release(fake(github=OSError('blocked')))
        self.assertEqual((r['version'], r['source']), ('0.7.10', 'website'))     # highest version of the website repo
        self.assertEqual(r['assets'], ['RCV-0.7.10-Linux-x86_64.rpm'])

    def test_offline(self):
        with self.assertRaises(RuntimeError) as cm:
            U.latest_release(fake(github=OSError('down'), stats=OSError('down')))
        self.assertIn('could not check for updates', str(cm.exception))


class TestSettingsAndAutoCheck(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        env = {'RCV_HOME': self.tmp.name}
        self.env = mock.patch.dict(os.environ, env)
        self.env.start()
        os.environ.pop('RCV_NO_UPDATE_CHECK', None)

    def tearDown(self):
        self.env.stop()
        self.tmp.cleanup()

    def test_default_off_from_source_and_on_when_installed(self):
        with mock.patch.object(U.P, 'frozen', return_value=False):
            self.assertFalse(U.settings()['auto'])
        with mock.patch.object(U.P, 'frozen', return_value=True):
            self.assertTrue(U.settings()['auto'])
            os.environ['RCV_NO_UPDATE_CHECK'] = '1'
            self.assertFalse(U.settings()['auto'])

    def test_auto_check_once_a_day_and_skip(self):
        U.save(auto=True)
        kw = dict(current='0.7.4', fetch=fake(), key='windows-x64', installed=True)
        r = U.auto_check(now=1e9, **kw)
        self.assertTrue(r['newer'])
        self.assertFalse(r['skipped'])
        self.assertIsNone(U.auto_check(now=1e9 + 3600, **kw))                 # checked an hour ago
        U.save(skip='0.8.0')
        r = U.auto_check(now=1e9 + U.AUTO_INTERVAL + 1, **kw)
        self.assertTrue(r['skipped'])                                           # the user skipped this version
        U.save(auto=False)
        self.assertIsNone(U.auto_check(now=1e12, **kw))
        saved = json.loads((Path(self.tmp.name) / 'settings.json').read_text())
        self.assertEqual((saved['update_check'], saved['update_skip']), (False, '0.8.0'))

    def test_auto_check_quiet_offline(self):
        U.save(auto=True)
        self.assertIsNone(U.auto_check(now=1e9, fetch=fake(github=OSError(), stats=OSError())))


class TestInterface(unittest.TestCase):
    def test_menu_dialog_and_translations(self):
        js = (ROOT / 'src/rcv/static/js/app.js').read_text(encoding='utf-8')
        self.assertIn("'Check for updates…'", js)
        self.assertIn('autoCheckForUpdates(this)', js)
        upd = (ROOT / 'src/rcv/static/js/update.js').read_text(encoding='utf-8')
        for code in ('zh-CN', 'de', 'ru'):
            s = json.loads((ROOT / f'src/rcv/static/i18n/{code}.json').read_text(encoding='utf-8'))['strings']
            for text in ('Check for updates…', 'Update available', 'Skip this version', 'Download',
                         'RCV {0} is available. You have version {1}.', 'Installer for this computer ({0}):'):
                self.assertIn(text, s, code)
        self.assertIn('Installer for this computer (', upd)


if __name__ == '__main__':
    unittest.main()
