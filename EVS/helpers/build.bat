x86_64-w64-mingw32-g++ -m64 -c libevs.cpp

x86_64-w64-mingw32-g++ -m64 -fPIC --shared -o libevs.so libevs.o