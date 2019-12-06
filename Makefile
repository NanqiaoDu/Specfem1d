cc = gcc

prom = specfem1d

deps = $(shell find include/*.h)
src = $(shell find src/*.c)

obj = $(src:%.c=%.o) 

$(prom): $(obj)
	$(cc) -o $(prom) $(obj)  -lm -lgsl -lgslcblas

%.o: %.c $(deps)
	$(cc) -g -c -Iinclude/ $< -o $@ 

clean:
	rm  $(obj)
