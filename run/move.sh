#!/bin/sh
DIRNAME="$1"
#rm -r signif/plot/$DIRNAME
mv plot/$DIRNAME signif/plot
mv outroot/$DIRNAME signif/outroot
