#!/bin/sh
mencoder "mf://$1/*.png" -mf fps=$3 -o $2 -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=9600
