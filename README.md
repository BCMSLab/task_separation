# task_separation

Code for the article "Cancer cells forgo translating the transcribed genes of non-specialized tasks"

To run this code,

- Build a docker image using the `Dockerfile` or pull the latest image [bcmslab/task_separation](https://hub.docker.com/r/bcmslab/task_separation)
- Run the three R scripts in `code/` in this order
    1. `code/prepare_data.R`: downloads the required data from source
    2. `code/transcription_translation.R`: excutes the code to perform the analysis and save the output
    3. `analysis.Rmd`: genrates the figures and graphs
