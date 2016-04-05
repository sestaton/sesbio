#!/bin/bash

cd `pwd`

jar=$HOME/apps/snpEff/snpEff.jar

java -Xmx4g -jar $jar eff \
     -noStats \
     -sequenceOntology \
     -hgvs GRCh37.75 \
     ../vcf2maf/data/test.vcf \
     > vcf2maf-test.snpeff.vcf
