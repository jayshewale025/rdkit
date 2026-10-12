"""Relocatable data paths for a wheel installation."""

import os

RDCodeDir = os.path.dirname(__file__)
_share = RDCodeDir
RDBaseDir = RDCodeDir
RDDataDir = os.path.join(_share, "Data")
RDDocsDir = os.path.join(_share, "Docs")
RDDemoDir = os.path.join(_share, "Demo")
RDBinDir = os.path.join(_share, "bin")
RDProjDir = os.path.join(_share, "Projects")
RDContribDir = os.path.join(_share, "Contrib")
