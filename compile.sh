
if [ ! -d $setupDir/build ]; then
  echo "build directory is not found, create ..."
  mkdir $setupDir/build
fi

cd build
if [ "$1" = "clean" ]; then
  echo "cleaning build directory ..."
  rm -r $setupDir/build/*
  echo "cmake -> $setupDir"
  cmake $setupDir
fi

if [ "$1" = "cmake" ]; then
  echo "cmake -> $setupDir"
  cmake $setupDir
fi

echo "make ..."
make
cd -
