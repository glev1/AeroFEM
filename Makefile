install: 
	pip install --upgrade pip &&\
		python3 -m pip install -r requirements.txt

lint:
	pylint --disable=R,C ./aerofem