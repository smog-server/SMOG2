#put your preferred java bin directory here
javabin=/Library/Java/JavaVirtualMachines/1.6.0.jdk/Contents/Commands
javashare=../../../share/java/

cp -r $javashare/argparser .
cp -r $javashare/ral .
cp -r $javashare/org/smogserver/io org/smogserver/io
cp -r $javashare/org/smogserver/util org/smogserver/util

$javabin/javac -target 1.6 -classpath . org/smogserver/scm/*.java org/smogserver/io/*.java org/smogserver/util/*.java org/smogserver/util/geom/*.java org/smogserver/util/math/*.java org/smogserver/util/exception/DoneException.java ral/*java

$javabin/jar vcmf manifest.txt SCM.jar  org/smogserver/scm/*.class argparser/*.class org/smogserver/io/*.class org/smogserver/util/*.class org/smogserver/util/geom/*.class org/smogserver/util/math/*.class org/smogserver/util/exception/DoneException.class

rm -r argparser
rm -r ral
rm -r org/smogserver/io
rm -r org/smogserver/util

mv SCM.jar ..

