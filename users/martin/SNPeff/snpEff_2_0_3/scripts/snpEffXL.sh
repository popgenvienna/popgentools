#!/bin/sh

DIR=$HOME/snpEff/

java -Xmx10G \
	-classpath "$DIR/lib/charts4j-1.2.jar:$DIR/lib/flanagan.jar:$DIR/lib/freemarker.jar:$DIR/lib/junit.jar:$DIR/lib/trove-2.1.0.jar:$DIR" \
	ca.mcgill.mcb.pcingola.snpEffect.commandLine.SnpEff \
	$*
