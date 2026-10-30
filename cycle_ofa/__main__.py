"""Support ``python -m cycle_ofa`` without installing the package."""

from .cli import main

raise SystemExit(main())
