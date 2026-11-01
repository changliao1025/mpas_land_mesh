#!/usr/bin/env python3
"""
Backward-compatibility driver for ocn30coast10lnd10.

Delegates directly to the config-driven workflow runner in run_workflow.py.
"""

from run_workflow import main

if __name__ == "__main__":
    main()
