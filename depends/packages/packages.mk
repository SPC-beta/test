packages:=boost openssl libevent gmp zlib backtrace tor bls-dash qrencode 
darwin_packages:=zeromq
linux_packages:=zeromq
native_packages :=

qt_linux_packages:=qt expat libxcb xcb_proto libXau xproto freetype fontconfig libxkbcommon libxcb_util libxcb_util_render libxcb_util_keysyms libxcb_util_image libxcb_util_wm

qrencode_linux_packages = qrencode
qrencode_android_packages = qrencode
qrencode_darwin_packages = qrencode
qrencode_mingw32_packages = qrencode

upnp_packages=miniupnpc

qt_darwin_packages=qt
qt_mingw32_packages=qt

bdb_packages=bdb

zmq_packages=zeromq

darwin_native_packages=
$(host_arch)_$(host_os)_native_packages+=native_b2

