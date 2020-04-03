if(`test -n "-L/home/lisa/miniconda3/lib -lxml2 -L/home/lisa/miniconda3/lib -lz -L/home/lisa/miniconda3/lib -llzma -lpthread -L/home/lisa/miniconda3/lib -lm -ldl"`) then

if(${?LD_LIBRARY_PATH}) then
    setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:-L/home/lisa/miniconda3/lib -lxml2 -L/home/lisa/miniconda3/lib -lz -L/home/lisa/miniconda3/lib -llzma -lpthread -L/home/lisa/miniconda3/lib -lm -ldl
else
   setenv LD_LIBRARY_PATH -L/home/lisa/miniconda3/lib -lxml2 -L/home/lisa/miniconda3/lib -lz -L/home/lisa/miniconda3/lib -llzma -lpthread -L/home/lisa/miniconda3/lib -lm -ldl
endif

endif
