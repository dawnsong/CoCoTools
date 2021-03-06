#!/usr/bin/env python
"""Script to commit the doc build outputs into the github-pages repo.

Use:

  gh-pages.py [tag]

If no tag is given, the current output of 'git describe' is used.  If given,
that is how the resulting directory will be named.

In practice, you should use either actual clean tags from a current build or
something like 'current' as a stable URL for the most current version of the """

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
from distutils.dir_util import copy_tree
import os
import re
import sys
from os import chdir as cd

from subprocess import Popen, PIPE, CalledProcessError, check_call

#-----------------------------------------------------------------------------
# Globals
#-----------------------------------------------------------------------------

pages_dir = 'website'
pages_branch = 'master'
html_dir = 'build/html'
pages_repo = 'git@github.com:cocotools/cocotools.github.io.git'

#-----------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------
def sh(cmd):
    """Execute command in a subshell, return status code."""
    return check_call(cmd, shell=True)


def sh2(cmd):
    """Execute command in a subshell, return stdout.

    Stderr is unbuffered from the subshell."""
    p = Popen(cmd, stdout=PIPE, shell=True)
    out = p.communicate()[0]
    retcode = p.returncode
    if retcode:
        raise CalledProcessError(retcode, cmd)
    else:
        return out.rstrip()


def sh3(cmd):
    """Execute command in a subshell, return stdout, stderr

    If anything appears in stderr, print it out to sys.stderr"""
    p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
    out, err = p.communicate()
    retcode = p.returncode
    if retcode:
        raise CalledProcessError(retcode, cmd)
    else:
        return out.rstrip(), err.rstrip()

def init_repo(path):
    """clone the gh-pages repo if we haven't already."""
    sh("git clone %s %s" % (pages_repo, path) )
    here = os.getcwdu()
    cd(path)
    sh('git checkout %s' % pages_branch)
    cd(here)

#-----------------------------------------------------------------------------
# Script starts
#-----------------------------------------------------------------------------
if __name__ == '__main__':
    # The tag can be given as a positional argument
    try:
        tag = sys.argv[1]
    except IndexError:
        try:
            tag = sh2('git describe --exact-match')
        except CalledProcessError:
            tag = "dev"   # Fallback

    startdir = os.getcwdu()
    if not os.path.exists(pages_dir):
        # init the repo
        init_repo(pages_dir)
    else:
        # ensure up-to-date before operating
        cd(pages_dir)
        sh('git checkout %s' % pages_branch)
        sh('git pull')
        cd(startdir)

    dest = pages_dir
    copy_tree(html_dir, dest)

    try:
        cd(pages_dir)
        stout = sh2('git status | head -1')
        branch = re.match('\# On branch (.*)$', stout).group(1)
        if branch != pages_branch:
            e = 'On %r, git branch is %r, MUST be "%s"' % \
                (pages_dir, branch, pages_branch)
            raise RuntimeError(e)

        sh('git add -A *')
        to_commit = sh2('git status -z')
        if not to_commit:
            print "\nNo new changes to commit and upload, nothing to do."
            sys.exit(0)
        else:
            sh('git commit -m"Updated doc release: %s"' % tag)
            print
            print 'Most recent 3 commits:'
            sys.stdout.flush()
            sh('git --no-pager log --oneline HEAD~3..')
    finally:
        cd(startdir)

    print
    print 'Now verify the build in: %r' % dest
    print "If everything looks good, 'git push'"
