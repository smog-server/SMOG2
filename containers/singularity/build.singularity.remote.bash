#!/bin/bash
# as long as the auth key is exported, then the git version of the singularity container can be updated with this:
tag=$1
scs-build build smog2.$tag.def library:smog-server/library/smog2:$tag
