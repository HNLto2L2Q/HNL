import FWCore.ParameterSet.Config as cms
import time

from functools import wraps
import os
import subprocess
import sys
import warnings


def diaper(f):
    # Decorator to make sure a function never crashes.
    @wraps(f)
    def wrapper(*args, **kwargs):
        try:
            return f(*args, **kwargs)
        except:
            # Just catch all the shit.
            print "version.py - exception caught and discarded"
            print "Unexpected error:", str(sys.exc_info()[0])
            return "UNKNOWN"
    return wrapper

_fwk_directory = os.path.join(
    os.environ['CMSSW_BASE'], 'src', 'HNL')

@diaper
def cmssw_version():
    return os.getenv('CMSSW_VERSION')

@diaper
def git_version():
    ''' Get the current commit hash of the framework, without using the git command

    Reads the information directly from the .git repository folder.
    '''
    HEAD_file = os.path.join(
        _fwk_directory, '.git', 'HEAD')
    if not os.path.exists(HEAD_file):
        warnings.warn("Could not extract git commit information!")
        return 'NO_IDEA'
    # Get current HEAD ref
    with open(HEAD_file, 'r') as head:
        head_ref = head.readline().split(':')[1].strip()
        commit_file = os.path.join(
            _fwk_directory, '.git', head_ref)
        # Read commit ID for current HEAD
        with open(commit_file, 'r') as commit:
            return commit.readline().strip()[0:7]

@diaper
def get_user():
    ''' Get the user name in a safe way '''
    if 'dboard_user' in os.environ:
        return os.environ['dboard_user']
    elif 'LOGNAME' in os.environ:
        return os.environ['LOGNAME']
    return 'UNKNOWN'

metaTree = cms.EDAnalyzer(
    'MetaNtuplizer',
    isMC=cms.bool(False),
    weightsSrc = cms.InputTag('externalLHEProducer'),
    commit=cms.string( git_version()),
    user=cms.string(   get_user()),
    cmsswVersion= cms.string( cmssw_version()),
    date=cms.string(time.strftime("%d %b %Y %H:%M:%S +0000", time.gmtime())),
    globalTag=cms.string(''),
    args = cms.string(''),
    hasLHE=cms.bool(True)
)
