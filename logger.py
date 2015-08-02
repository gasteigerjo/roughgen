from time import strftime

def log(message):
    timestamp = strftime("%Y-%m-%d %H:%M:%S")
    print(timestamp + ": " + message)
