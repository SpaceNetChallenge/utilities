"""Run pytest from pythonw for matplotlib backend errors"""
if __name__ == '__main__':
    import sys
    from py.test import main

    sys.exit(main())