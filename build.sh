~/bin/google-compiler --js=../starjs/lib/StarJs{.js,/Math.js,/Vector.js,/Time.js,/Coord.js,/Kepler.js,/Solar.js} --js=starmap.js --js=../starjs/lib/export.js --js_output_file=starmap.min.js --compilation_level ADVANCED_OPTIMIZATIONS
cat ../starjs/lib/StarJs{.js,/Math.js,/Vector.js,/Time.js,/Coord.js,/Kepler.js,/Solar.js} starmap.js > starmap.debug.js
