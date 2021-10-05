#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
import copy

from .Fastb import Fastb

BASE="/Users/bmajoros/python/test/data"
#filename=BASE+"/iter0_peak1858.standardized_across_all_timepoints.t05.fastb"
filename=BASE+"/test.fastb"

fastb=Fastb(filename)
fastb.renameTrack("EP300","p300")
p300=fastb.getTrackByName("p300")
p300_2=copy.deepcopy(p300)
p300_2.id="newtrack"
fastb.addTrack(p300_2)
fastb.dropTrack("DNase")

numTracks=fastb.numTracks()
print(numTracks,"tracks")
for i in range(numTracks):
    track=fastb.getIthTrack(i)
    print(track.getID(),track.getLength())

fastb=fastb.slice(500,1500)
fastb.save(BASE+"/test-out.fastb")
