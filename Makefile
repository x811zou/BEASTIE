.PHONY: clean install

PYTHON_FILES=$(wildcard BEASTIE/*/*.py)

clean:
	-rm -rf BEASTIE.egg-info build dist

docker: clean $(PYTHON_FILES)
	docker build .

dist: $(PYTHON_FILES)
	if [ -z $$(pip list | grep -e "^build\s") ]; then pip install build; fi
	python -m build

install: dist
	pip install dist/*.whl
