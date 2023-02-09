#!/bin/bash
# as long as the auth key is exported, then the git version of the singularity container can be updated with this:
tag=$1
scs-build build smogopensmog.$tag.def library:smog-server/library/smogopensmog:$tag
