#put your preferred java bin directory here
javabin=/Library/Java/JavaVirtualMachines/1.6.0.jdk/Contents/Commands
javashare=../../../share/java/

cp -r $javashare/argparser .
cp -r $javashare/ral .
cp -r $javashare/org/smogserver/io org/smogserver/
cp -r $javashare/org/smogserver/util org/smogserver/

$javabin/javac -target 1.6 -classpath . org/smogserver/wham/*.java org/smogserver/io/*.java org/smogserver/util/*.java org/smogserver/util/geom/*.java org/smogserver/util/math/*.java org/smogserver/util/exception/DoneException.java ral/*java &> error

if [ $( grep error error | wc | awk '{print $1}' ) -gt 0 ]
then
echo there were compilation errors
echo -----------------------------
cat error
exit
fi 

$javabin/jar vcmf manifest.txt WHAM.jar org/smogserver/wham/*.class argparser/*.class org/smogserver/io/*.class org/smogserver/util/*.class org/smogserver/util/geom/*.class org/smogserver/util/math/*.class org/smogserver/util/exception/DoneException.class ral/*class

rm -r argparser
rm -r ral
rm -r org/smogserver/io
rm -r org/smogserver/util

mv WHAM.jar ..
rm error
