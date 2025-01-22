import os


def isfileempty(f):
    """Checks if the opened file is empty. To be used in file context manager"""
    f.seek(0, 2)  # Move to the end of the file
    if f.tell() > 0:  # If the file is not empty
        return False
    else:
        return True


def get_base_directory(base_folder_name, current_path):
    try:
        while True:
            # Check if the marker file/directory exists
            if os.path.exists(os.path.join(current_path, base_folder_name)):
                break
                # Move up one directory
            parent_path = os.path.dirname(current_path)
            if parent_path == current_path:  # We've reached the root directory
                raise Exception()
            current_path = parent_path

    except Exception:
        print(f"{base_folder_name} not found in any parent directory.")

    finally:
        return os.path.join(current_path, base_folder_name)


def is_not_ragged(array):
    # Check if the array is empty
    if not array:
        return True  # An empty array is not considered ragged

    # Get the length of the first inner list
    first_length = len(array[0])

    # Check if all inner lists are of the same length
    for inner_list in array:
        if len(inner_list) != first_length:
            return False  # Found an inner list with a different length

    return True  # All inner lists have the same length

def single_retry(func):
    """Wrapper to single-retry the given function."""
    def wrapper(*args, **kwargs):
        try:
            # Try to call the function
            return func(*args, **kwargs)
        except AssertionError:
            print("First attempt failed. Retrying once...")
            # Try again
            try:
                return func(*args, **kwargs)
            except AssertionError:
                print("Second attempt failed with error. No more retries.")
                return func(*args, **kwargs)
    return wrapper