all: larva

larva: larva.cpp p.value.calc.cpp moment.estimator.cpp
	./dependency-check.sh && g++ -Wall -o larva larva.cpp p.value.calc.cpp moment.estimator.cpp
