#!/usr/bin/env python3
"""
ImmunoAtlas — gzip-aware HTTP server
Usage: python3 server.py [port]
"""
import os, sys, gzip
from http.server import HTTPServer, SimpleHTTPRequestHandler

PORT = int(sys.argv[1]) if len(sys.argv) > 1 else 8080
ROOT = os.path.dirname(os.path.abspath(__file__))

class GzipHandler(SimpleHTTPRequestHandler):
    """Serve pre-compressed .gz files when the browser accepts gzip."""

    def do_GET(self):
        accept_enc = self.headers.get('Accept-Encoding', '')
        if 'gzip' in accept_enc:
            # strip query string for file lookup
            path = self.path.split('?')[0].split('#')[0]
            fs_path = os.path.join(ROOT, path.lstrip('/'))
            gz_path = fs_path + '.gz'
            if os.path.isfile(gz_path):
                self._serve_gz(gz_path, fs_path)
                return
        super().do_GET()

    def _serve_gz(self, gz_path, orig_path):
        """Send the pre-gzipped file with appropriate headers."""
        stat = os.stat(gz_path)
        size = stat.st_size

        # Guess content type from the *original* filename
        import mimetypes
        ctype, _ = mimetypes.guess_type(orig_path)
        if ctype is None:
            ctype = 'application/octet-stream'

        self.send_response(200)
        self.send_header('Content-Type', ctype)
        self.send_header('Content-Encoding', 'gzip')
        self.send_header('Content-Length', str(size))
        self.send_header('Cache-Control', 'no-cache')
        self.end_headers()

        with open(gz_path, 'rb') as f:
            while True:
                chunk = f.read(65536)
                if not chunk:
                    break
                try:
                    self.wfile.write(chunk)
                except BrokenPipeError:
                    break

    def log_message(self, fmt, *args):
        # Show file size info for large files
        super().log_message(fmt, *args)


if __name__ == '__main__':
    os.chdir(ROOT)
    server = HTTPServer(('0.0.0.0', PORT), GzipHandler)
    print(f'ImmunoAtlas server running at http://0.0.0.0:{PORT}/')
    print(f'Root: {ROOT}')
    print('Press Ctrl+C to stop.')
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print('\nStopped.')
