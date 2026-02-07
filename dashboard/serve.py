#!/usr/bin/env python3
"""
Simple HTTP server for Pycnopodia dashboard.
Run: python3 serve.py
Then open: http://localhost:8090
"""

import http.server
import socketserver
import os
import sys

PORT = 8090

class MyHTTPRequestHandler(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        # Enable CORS for local development
        self.send_header('Access-Control-Allow-Origin', '*')
        super().end_headers()

def main():
    # Change to dashboard directory
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    
    # Check if data.json exists
    if not os.path.exists('data.json'):
        print("⚠️  data.json not found!")
        print("Please run: python3 ../run_pacific_coast.py --json-only")
        sys.exit(1)
    
    Handler = MyHTTPRequestHandler
    
    with socketserver.TCPServer(("", PORT), Handler) as httpd:
        print(f"🌊 Pycnopodia Dashboard Server")
        print(f"━" * 50)
        print(f"Serving at: http://localhost:{PORT}")
        print(f"Press Ctrl+C to stop")
        print(f"━" * 50)
        
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print("\n\n👋 Server stopped")
            sys.exit(0)

if __name__ == "__main__":
    main()
