"""Entry point of the RCV desktop builds (PyInstaller).

The same executable also runs the prediction engine (``RCV engine predict ...``) and its worker
processes, so ``freeze_support`` must run before anything else.
"""
import multiprocessing
import sys

if __name__ == '__main__':
    multiprocessing.freeze_support()
    from rcv.app import main
    sys.exit(main())
