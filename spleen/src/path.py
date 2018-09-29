# -*- coding: utf-8 -*-

u""" path.py - An object representing a path to a file or directory.

Example:

from path import path
d = Path('/home/guido/bin')
for f in d.files('*.py'):
    f.chmod(0755)

This module requires Python 2.2 or later.


URL:     http://www.jorendorff.com/articles/python/path
Author:  Jason Orendorff <jason@jorendorff.com> (and others - see the url!)
Date:    7 Mar 2004

Adapted for stdlib by: Reinhold Birkenfeld, July 2005
Modified by Björn Lindqvist <bjourne@gmail.com>, January 2006
"""

__all__=[
           u'DirPath',
           u'FilePath',
           u'Path',
           u'format_path'
           ]
def tr_(chaine):return chaine
try : from PyQt4.QtCore import QString
except ImportError : QString=str
# TODO
#   - Better error message in listdir() when self isn't a
#     directory. (On Windows, the error message really sucks.)
#   - Make sure everything has a good docstring.
#   - Add methods for regex find and replace.
#   - Perhaps support arguments to touch().
#   - Could add split() and join() methods that generate warnings.
#   - Note:  __add__() technically has a bug, I think, where
#     it doesn't play nice with other types that implement
#     __radd__().  Test this.

import fnmatch
import glob
import os
import shutil
import locale
import string

#from _concha.core.system import error, uprint, ENCODING

__version__=u'2.0.4'
ENCODING='utf-8'
# Universal newline support
_textmode=u'r'
if hasattr(file,u'newlines'):
    _textmode=u'U'

_base=unicode
GOOD_TYPES=(unicode,QString,str)
def goodInitType(arg):
    r=False
    for t in GOOD_TYPES :
        r=r or isinstance(arg,t)
        if r : return True
    return False

class Path(_base):
    u""" Represents a filesystem path.

    For documentation on individual methods, consult their
    counterparts in os.path.
    """
    special_chars=u'"!'
    # --- Special Python methods.
    def __new__(cls,*args):
        u"""
        Creates a new path object concatenating the *args.  *args
        may only contain Path objects or strings.  If *args is
        empty, an exception is raised.
        """
        if not args:
            pass
#            raise ValueError, u'At least one argument required'
        else :
            pass

        for arg in args:
#            if not isinstance(arg, unicode):
            if not goodInitType(arg) :
                raise TypeError,tr_(u'E006 - Path, unicode expected, got %(TYPE)s'%dict(
                  TYPE=type(arg).__name__))

        if len(args)<=1:
#          if not args[0].strip() :
#            raise ValueError, tr_(u'E007 - Dangerous value, you must use u"." for current directory')
            return unicode.__new__(cls,*args)
        return cls(os.path.join(*args))

    def __repr__(self):
        return u'%s(%r)'%(self.__class__.__name__,unicode(self))

    def __str__ (self):
        return self.encode(ENCODING)

    # Adding a path and a string yields a path.
    def __add__(self,more):
        return self.__class__(unicode.__add__(unicode(self),unicode(more)))

    def __radd__(self,other):
        return self.__class__(other+_base(self))

    @classmethod
    def cwd(cls):
        u""" Returns the current working directory as a path object. """
        return Path(os.getcwdu())

    # --- Operations on path strings.
    def join(self,seq):
        raise NotImplementedError,u"To join path : Path(p1, p2, p3, ...)"

    def _abspath(self):
        # see: http://bugs.python.org/issue3426
        path=os.path.abspath(self.encode(ENCODING)).decode(ENCODING)
        return self.__class__(path)

    def shellpath(self):
        u"""
        Returns a path directly usable in a shell.
        /home/guillaume/mon "dossier" tordu
        -> "/home/guillaume/mon\ \"dossier\"\ tordu
        """
        path=unicode(self)
        for c in list(self.special_chars) :
            path=path.replace(c,u'\%s'%c)
        return '"%s"'%path

    def normcase(self):
        return self.__class__(os.path.normcase(self))

    def normpath(self):
        # see: http://bugs.python.org/issue3426
        path=os.path.normpath(self.encode(ENCODING)).decode(ENCODING)
        return self.__class__(path)

    def realpath(self):
        # see: http://bugs.python.org/issue3426
        path=os.path.realpath(self)
        if isinstance(path,str) :
            path=path.decode(ENCODING)
        return self.__class__(path)

    def splitext(self):
        return os.path.splitext(self)

    def replaceext(self,ext,old_ext=None):
        u"""
        Changes current extension to ext.
        If extension contains more than one dot (ex: .tar.gz) you can specify
        it with "old_ext" argument.
        
        >>> p = Path(u'path.py')
        >>> p.replaceext(u'.rst')
        Path(u'path.rst')
        
        >>> p = Path(u'concha.tar.gz')
        >>> p.replaceext(u'.tgz')
        Path(u'concha.tar.tgz')
        >>> p.replaceext(u'.tgz', u'.tar.gz')
        Path(u'concha.tgz')
        """
        if old_ext is None :
            path,old_ext=self.splitext()
            return self.__class__(path+ext)
        else :
            path_ext=self[-len(old_ext):]
            if path_ext!=old_ext :
                raise ValueError,'old_ext %(OLD_EXT)r and path ext %(PATH_EXT)r do not match.'%dict(
                                   OLD_EXT=old_ext,PATH_EXT=path_ext)
            else :
                return self.__class__(self[:-len(old_ext)]+ext)

    def expanduser(self):
        return self.__class__(os.path.expanduser(self))

    def expandvars(self):
        return self.__class__(os.path.expandvars(self))

    def expand(self):
        u""" Cleans up a filename by calling expandvars(),
        expanduser(), and normpath() on it.

        This is commonly everything needed to clean up a filename
        read from a configuration file, for example.
        """
        return self.expandvars().expanduser().normpath()

    def _get_namebase(self):
        base,ext=os.path.splitext(self.name)
        return base

    def _get_ext(self):
        f,ext=os.path.splitext(_base(self))
        return ext

    def _get_drive(self):
        drive,r=os.path.splitdrive(self)
        return self.__class__(drive)

    def _get_dirname(self):
        return self.__class__(os.path.dirname(self))

    parent=property(
        _get_dirname,None,None,
        u""" This path's parent directory, as a new path object.

        For example, Path('/usr/local/lib/libpython.so').parent == Path('/usr/local/lib')
        """)
    dirname=parent

    name=property(
        os.path.basename,None,None,
        u""" The name of this file or directory without the full path.

        For example, Path('/usr/local/lib/libpython.so').name == 'libpython.so'
        """)

    namebase=property(
        _get_namebase,None,None,
        u""" The same as path.name, but with one file extension stripped off.

        For example, Path('/home/guido/python.tar.gz').name     == 'python.tar.gz',
        but          Path('/home/guido/python.tar.gz').namebase == 'python.tar'
        """)

    ext=property(
        _get_ext,None,None,
        u""" The file extension, for example '.py'. """)

    drive=property(
        _get_drive,None,None,
        u""" The drive specifier, for example 'C:'.
        This is always empty on systems that don't use drive specifiers.
        """)
    abspath=property(_abspath,None,None,u"The absolute path")

    def splitpath(self):
        u""" p.splitpath() -> Return (p.parent, p.name). """
        parent,child=os.path.split(self)
        return self.__class__(parent),child

    def stripext(self):
        u""" p.stripext() -> Remove one file extension from the path.

        For example, Path('/home/guido/python.tar.gz').stripext()
        returns Path('/home/guido/python.tar').
        """
        return Path(os.path.splitext(self)[0])

    if hasattr(os.path,u'splitunc'):
        def splitunc(self):
            unc,rest=os.path.splitunc(self)
            return self.__class__(unc),rest

        def _get_uncshare(self):
            unc,r=os.path.splitunc(self)
            return self.__class__(unc)

        uncshare=property(
            _get_uncshare,None,None,
            u""" The UNC mount point for this path.
            This is empty for paths on local drives. """)

    def splitall(self):
        u""" Returns a list of the path components in this path.

        The first item in the list will be a path.  Its value will be
        either os.curdir, os.pardir, empty, or the root directory of
        this path (for example, '/' or 'C:\\').  The other items in
        the list will be strings.

        path.path(*result) will yield the original path.
        """
        parts=[]
        loc=self
        while loc!=os.curdir and loc!=os.pardir:
            prev=loc
            loc,child=prev.splitpath()
            loc=self.__class__(loc)
            if loc==prev:
                break
            parts.append(child)
        parts.append(loc)
        parts.reverse()
        return parts

    def shortestpath(self):
        u"""
        Returns the shortest path between _abspath and relpath.
        This form is useful for display. 
        """
        relpath=len(unicode(self.relpath()))
        _abspath=len(unicode(self._abspath()))
        if (_abspath<=relpath) :
            return self._abspath()
        else :
            return self.relpath()

    def relpath(self):
        u""" Returns this path as a relative path,
        based from the current working directory.
        """
        return self.__class__.cwd().relpathto(self)

    def relpathto(self,dest):
        u""" Returns a relative path from self to dest.

        If there is no relative path from self to dest, for example if
        they reside on different drives in Windows, then this returns
        dest._abspath().
        """
        if not isinstance(dest,unicode):
            raise TypeError,tr_(u'E006 - Path, unicode expected, got %(TYPE)s'%dict(
                TYPE=type(dest).__name__))

        origin=self._abspath()
        dest=self.__class__(dest)._abspath()

        orig_list=origin.normcase().splitall()
        # Don't normcase dest!  We want to preserve the case.
        dest_list=dest.splitall()

        if orig_list[0]!=os.path.normcase(dest_list[0]):
            # Can't get here from there.
            return dest

        # Find the location where the two paths start to differ.
        i=0
        for start_seg,dest_seg in zip(orig_list,dest_list):
            if start_seg!=os.path.normcase(dest_seg):
                break
            i+=1

        # Now i is the point where the two paths diverge.
        # Need a certain number of "os.pardir"s to work up
        # from the origin to the point of divergence.
        segments=[os.pardir.decode(ENCODING)]*(len(orig_list)-i)
        # Need to add the diverging part of dest_list.
        segments+=dest_list[i:]
        if len(segments)==0:
            # If they happen to be identical, use os.curdir.
            return self.__class__(os.curdir.decode(ENCODING))
        else:
            return self.__class__(os.path.join(*segments))


    # --- Listing, searching, walking, and matching

    def listdir(self,pattern=None):
        u""" D.listdir() -> Lists of items in this directory.

        Use D.files() or D.dirs() instead if you want a listing
        of just files or just subdirectories.

        The elements of the list are path objects.

        With the optional 'pattern' argument, this only lists
        items whose names match the given pattern.
        """
        names=os.listdir(unicode(self))
        if pattern is not None:
            names=fnmatch.filter(names,pattern)
        return [Path(unicode(self),child) for child in names]

    def dirs(self,pattern=None):
        u""" D.dirs() -> List of this directory's subdirectories.

        The elements of the list are path objects.
        This does not walk recursively into subdirectories
        (but see path.walkdirs).

        With the optional 'pattern' argument, this only lists
        directories whose names match the given pattern.  For
        example, d.dirs('build-*').
        """
        return [p for p in self.listdir(pattern) if p.isdir()]

    def files(self,pattern=None):
        u""" D.files() -> List of the files in this directory.

        The elements of the list are path objects.
        This does not walk into subdirectories (see path.walkfiles).

        With the optional 'pattern' argument, this only lists files
        whose names match the given pattern.  For example,
        d.files('*.pyc').
        """

        return [p for p in self.listdir(pattern) if p.isfile()]

    def walk(self,pattern=None):
        u""" D.walk() -> iterator over files and subdirs, recursively.

        The iterator yields path objects naming each child item of
        this directory and its descendants.  This requires that
        D.isdir().

        This performs a depth-first traversal of the directory tree.
        Each directory is returned just before all its children.
        """
        for child in self.listdir():
            if pattern is None or child.match(pattern):
                yield child
            if child.isdir():
                for item in child.walk(pattern):
                    yield item

    def walkdirs(self,pattern=None):
        u""" D.walkdirs() -> iterator over subdirs, recursively.

        With the optional 'pattern' argument, this yields only
        directories whose names match the given pattern.  For
        example, mydir.walkdirs('*test') yields only directories
        with names ending in 'test'.
        """
        for child in self.dirs():
            if pattern is None or child.match(pattern):
                yield child
            for subsubdir in child.walkdirs(pattern):
                yield subsubdir

    def walkfiles(self,pattern=None):
        u""" D.walkfiles() -> iterator over files in D, recursively.

        The optional argument, pattern, limits the results to files
        with names that match the pattern.  For example,
        mydir.walkfiles('*.tmp') yields only files with the .tmp
        extension.
        """
        for child in self.listdir():
            #if child.isfile():
            if os.path.isfile(child):
                if pattern is None or child.match(pattern):
                    yield child
            #elif child.isdir():
            elif os.path.isdir(child):
                for f in child.walkfiles(pattern):
                    yield f

    def match(self,pattern):
        u""" Return True if self.name matches the given pattern.

        pattern - A filename pattern with wildcards,
        for example '*.py'.
        """
        return fnmatch.fnmatch(self.name,pattern)

    def matchcase(self,pattern):
        u""" Test whether the path matches pattern, returning true or
        false; the comparison is always case-sensitive.
        """
        return fnmatch.fnmatchcase(self.name,pattern)

    def glob(self,pattern):
        u""" Return a list of path objects that match the pattern.

        pattern - a path relative to this directory, with wildcards.

        For example, Path('/users').glob('*/bin/*') returns a list
        of all the files users have in their bin directories.
        """
        return map(Path,glob.glob(_base(Path(self,pattern))))

    # --- Methods for querying the filesystem.

    exists=os.path.exists
    isabs=os.path.isabs
    isdir=os.path.isdir
    isfile=os.path.isfile
    islink=os.path.islink
    ismount=os.path.ismount

    if hasattr(os.path,u'samefile'):
        samefile=os.path.samefile

    def atime(self):
        u"""Last access time of the file."""
        return os.path.getatime(self)

    def mtime(self):
        u"""Last-modified time of the file."""
        return os.path.getmtime(self)

    def ctime(self):
        u"""
        Return the system's ctime which, on some systems (like Unix)
        is the time of the last change, and, on others (like Windows),
        is the creation time for path.

        The return value is a number giving the number of seconds
        since the epoch (see the time module). Raise os.error if the
        file does not exist or is inaccessible.
        """
        return os.path.getctime(self)

    def size(self):
        u"""Size of the file, in bytes."""
        return os.path.getsize(self)

    if hasattr(os,u'access'):
        def access(self,mode):
            u""" Return true if current user has access to this path.

            mode - One of the constants os.F_OK, os.R_OK, os.W_OK, os.X_OK
            """
            return os.access(self,mode)

    def stat(self):
        u""" Perform a stat() system call on this path. """
        return os.stat(self)

    def lstat(self):
        u""" Like path.stat(), but do not follow symbolic links. """
        return os.lstat(self)

    if hasattr(os,u'statvfs'):
        def statvfs(self):
            u""" Perform a statvfs() system call on this path. """
            return os.statvfs(self)

    if hasattr(os,u'pathconf'):
        def pathconf(self,name):
            return os.pathconf(self,name)


    # --- Modifying operations on files and directories

    def utime(self,times):
        u""" Set the access and modified times of this file. """
        os.utime(self,times)

    def chmod(self,mode):
        os.chmod(self,mode)

    if hasattr(os,u'chown'):
        def chown(self,uid,gid):
            os.chown(self,uid,gid)

    def rename(self,new):
        os.rename(self,new)

    def renames(self,new):
        os.renames(self,new)


    # --- Create/delete operations on directories

    def mkdir(self,mode=0777):
        os.mkdir(self,mode)

    def makedirs(self,mode=0777):
        os.makedirs(self,mode)

    def rmdir(self):
        os.rmdir(self)

    def removedirs(self):
        os.removedirs(self)


    # --- Modifying operations on files

    def touch(self):
        u""" Set the access/modified times of this file to the current time.
        Create the file if it does not exist.
        """
        fd=os.open(self,os.O_WRONLY|os.O_CREAT,0666)
        os.close(fd)
        os.utime(self,None)

    def remove(self):
        os.remove(self)

    def unlink(self):
        os.unlink(self)


    # --- Links

    if hasattr(os,u'link'):
        def link(self,newpath):
            u""" Create a hard link at 'newpath', pointing to this file. """
            os.link(self,newpath)

    if hasattr(os,u'symlink'):
        def symlink(self,newlink):
            u""" Create a symbolic link at 'newlink', pointing here. """
            os.symlink(self,newlink)

    if hasattr(os,u'readlink'):
        def readlink(self):
            u""" Return the path to which this symbolic link points.

            The result may be an absolute or a relative path.
            """
            return self.__class__(os.readlink(self))

        def readlinkabs(self):
            u""" Return the path to which this symbolic link points.

            The result is always an absolute path.
            """
            p=self.readlink()
            if p.isabs():
                return p
            else:
                return self.__class__(self.parent,p)._abspath()

    # --- High-level functions from shutil

    copyfile=shutil.copyfile
    copymode=shutil.copymode
    copystat=shutil.copystat
    copy=shutil.copy
    copy2=shutil.copy2
    copytree=shutil.copytree
    if hasattr(shutil,u'move'):
        move=shutil.move
    rmtree=shutil.rmtree

    # --- Special stuff from os

    if hasattr(os,u'chroot'):
        def chroot(self):
            os.chroot(self)

class FilePath (Path):
    u"""
    Restricted Path for files only
    """

class DirPath (Path):
    u"""
    Restricted Path for directories only
    """

def format_path (path,_abspath=True,expanduser=True,pathbook=None,pathid=None):
    if path is None or path is u'' :
        if pathid is None or pathbook is None:
            return None
        else :
            path=pathbook(pathid)
    path=Path(path)
    if expanduser : path=path.expanduser()
    if _abspath : path=path._abspath()
    return path

if __name__=="__main__":
    p=Path('.')
    print p,p.abspath
