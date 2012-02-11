#!/usr/bin/env python
#-*- coding: UTF-8 -*-

import urllib
import os

class URLOpener(urllib.FancyURLopener):
    def http_error_206(self, url, fp, errcode, errmsg, headers, data=None):
        pass

def downloader(url, dest, retries=5):
    if retries == 0:
        return False
    opener = URLOpener()
    connection = opener.open(url)
    totalsize = int(connection.headers['Content-Length'])

    if os.path.exists(dest):
        output = open(dest, "ab")
        currentsize = os.path.getsize(dest)
        opener.addheader("Range","bytes=%s-" % currentsize)
        connection = opener.open(url)
    else:
        output = open(dest, "wb")
        currentsize = 0

    while currentsize < totalsize:
        data = connection.read(81920)
        output.write(data)
        currentsize += len(data)
        print "%d, " % (100 * currentsize / totalsize),
    
    if totalsize == currentsize:
        print ""
        return True
    else:
        return downloader(url, dest, retries - 1)
