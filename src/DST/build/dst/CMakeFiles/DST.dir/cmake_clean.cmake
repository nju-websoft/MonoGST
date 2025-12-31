file(REMOVE_RECURSE
  "libDST.a"
  "libDST.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/DST.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
