!#/bin/bash
buildDirectory='./build'

#Cleanup the build directory
if [ -d "$buildDirectory" ]; then
  rm -rf $buildDirectory
fi

#create a new build directory
mkdir $buildDirectory

#enter the build directory
cd $buildDirectory

# launch cmake
cmake ..

# launch make
make



