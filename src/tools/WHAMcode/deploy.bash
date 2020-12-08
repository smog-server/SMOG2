#put your preferred java bin directory here
javabin=/Library/Java/JavaVirtualMachines/1.6.0.jdk/Contents/Commands

$javabin/javac -target 1.6 -classpath . noel/folding/wham/*.java noel/io/*.java noel/util/*.java noel/util/geom/*.java noel/util/math/*.java noel/util/exception/DoneException.java ral/*java

$javabin/jar vcmf manifest.txt WHAM.jar  noel/folding/wham/*.class argparser/*.class noel/io/*.class noel/util/*.class noel/util/geom/*.class noel/util/math/*.class noel/util/exception/DoneException.class ral/*class

cp WHAM.jar ..

