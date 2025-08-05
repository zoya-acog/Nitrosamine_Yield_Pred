build-image:
	docker build -t nitrosamine-yield-pred .
run-container:
	docker run --gpus all -d --label description="bhawakshi's nitrosamine container" --name nitrosamine -v $(PWD):/home/jupyter/ nitrosamine-yield-pred

stop-container:
	docker stop nitrosamine
remove-container:
	docker rm -f nitrosamine
