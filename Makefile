.PHONY: setup pipeline dashboard

setup:
	pip install -r requirements.txt

pipeline:
	python load_data.py
	python run_pipeline.py

dashboard:
	streamlit run dashboard.py
