.PHONY: clean install

PYTHON_FILES=$(wildcard BEASTIE/*/*.py)

clean:
	-rm -rf BEASTIE.egg-info build dist

docker: clean $(PYTHON_FILES)
	docker build .

dist: $(PYTHON_FILES)
	if [ -z $$(pip3 list | grep -e "^build\s") ]; then pip3 install build; fi
	python3 -m build

install: dist
	pip3 install dist/*.whl
