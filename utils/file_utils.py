import os

def set_tests_dir():

    # Check if we are already in 'Tests' directory
    if os.path.basename(os.getcwd()) == "Tests":
        print("Already in the 'Tests' directory!")
        print(f"Current working directory: {os.getcwd()}")
    else:
        # Check if 'Tests' exists in current directory
        if not os.path.exists("Tests"):
            os.makedirs("Tests")
        os.chdir("Tests")
        print(f"Moved into 'Tests' directory. Current working directory: {os.getcwd()}")