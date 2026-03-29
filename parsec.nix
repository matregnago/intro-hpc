{ stdenv
, fetchgit
, cmake
, pkg-config
, hwloc
, openmpi
, gfortran
, flex
, bison
, python3
}:

stdenv.mkDerivation rec {
  pname = "parsec";
  version = "mymaster";

src = fetchgit {
  url = "https://bitbucket.org/mfaverge/parsec.git";
  rev = "6022a61dc96c25f11dd2aeabff2a5b3d7bce867d";
  fetchSubmodules = true;
  hash = "sha256-iJEIJLdakIxGFz+qXNTcbGOISzXbI/XLVq45k4GqhoM=";
};

  nativeBuildInputs = [
    cmake
    pkg-config
    gfortran
    flex
    bison
    python3
  ];

  buildInputs = [
    hwloc
    openmpi
    stdenv.cc.cc.lib
  ];

cmakeFlags = [
  "-DBUILD_SHARED_LIBS=ON"
  "-DCMAKE_POLICY_VERSION_MINIMUM=3.5"
];

env.NIX_CFLAGS_COMPILE = toString [
  "-Wno-error=implicit-function-declaration"
  "-Wno-error=incompatible-pointer-types"
  "-Wno-implicit-function-declaration"
  "-Wno-incompatible-pointer-types"
];
}