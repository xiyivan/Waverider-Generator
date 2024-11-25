import tkinter as tk
from tkinter import ttk
import numpy as np
import subprocess
import threading
import os
import sys

# Function to run the user_input.py script
def run_script():
    # Check and change the working directory to the current directory
    current_dir = os.path.dirname(os.path.abspath(__file__))
    if os.getcwd() != current_dir:
        os.chdir(current_dir)
    
    # Create the temp directory if it doesn't exist
    if not os.path.exists("temp"):
        os.makedirs("temp")
    
    with open("user_input.py", "w") as f:
        f.write(f"""
import numpy as np

#   Basic Flowfield
M1 = {M1_var.get()}
gamma = {gamma_var.get()}
beta = {beta_var.get()}
beta_rad = np.radians(beta)

#   Define TE constants
L = {L_var.get()}                                # Max. Length (See User Menu)
Rs = L * np.tan(np.radians(beta))    # Trailing Edge Parameter (See User Menu)
R1 = 0.2 * Rs                        # Trailing Edge Parameter (See User Menu)
W2 = 0.8 * Rs                        # Trailing Edge Parameter (See User Menu)

#   Resolution
N = {N_var.get()}                      # Leading Edge Resolution (500 or more is recommended!)
N_l = {N_l_var.get()}                     # Lower Surface Resolution
N_up = {N_up_var.get()}                    # Upper Surface Resolution

#   Trailing Edge Curve
a = -R1
b = 2 * (R1 - np.sqrt(Rs**2 - W2**2)) / W2**2
c = (np.sqrt(Rs**2 - W2**2) - R1) / W2**4

#   Trailing Edge Function
def z(y):
    func = {trailing_edge_function_var.get()}    # Define Your Own Trailing Edge Function
    return func
""")
    # Redirect stdout and stderr to the text widget
    sys.stdout = TextRedirector(output_text, "stdout")
    sys.stderr = TextRedirector(output_text, "stderr")
    
    subprocess.run(["python", "user_input.py"])
    subprocess.run(["python", "main.py"])
    update_image()

# Function to update the image
def update_image():
    img_path = "./temp/cone_fig.png"
    if os.path.exists(img_path):
        img = tk.PhotoImage(file=img_path)
        image_label.config(image=img)
        image_label.image = img

# Class to redirect stdout and stderr to the text widget
class TextRedirector(object):
    def __init__(self, widget, tag="stdout"):
        self.widget = widget
        self.tag = tag

    def write(self, str):
        self.widget.insert(tk.END, str, (self.tag,))
        self.widget.see(tk.END)

    def flush(self):
        pass

# Create the main window
root = tk.Tk()
root.title("GUI Program")

# Create the user input frame
input_frame = ttk.Frame(root)
input_frame.grid(row=0, column=0, padx=10, pady=10, sticky="nw")

# Basic Flowfield
ttk.Label(input_frame, text="Basic Flowfield").grid(row=0, column=0, columnspan=2)
M1_var = tk.DoubleVar(value=8)
ttk.Label(input_frame, text="M1").grid(row=1, column=0)
ttk.Entry(input_frame, textvariable=M1_var).grid(row=1, column=1)
gamma_var = tk.DoubleVar(value=1.4)
ttk.Label(input_frame, text="gamma").grid(row=2, column=0)
ttk.Entry(input_frame, textvariable=gamma_var).grid(row=2, column=1)
beta_var = tk.DoubleVar(value=16.5)
ttk.Label(input_frame, text="beta").grid(row=3, column=0)
ttk.Entry(input_frame, textvariable=beta_var).grid(row=3, column=1)

# Define TE constants
ttk.Label(input_frame, text="Define TE constants").grid(row=4, column=0, columnspan=2)
L_var = tk.DoubleVar(value=2)
ttk.Label(input_frame, text="L").grid(row=5, column=0)
ttk.Entry(input_frame, textvariable=L_var).grid(row=5, column=1)
Rs_var = tk.DoubleVar(value=2 * np.tan(np.radians(16.5)))
R1_var = tk.DoubleVar(value=0.2 * Rs_var.get())
W2_var = tk.DoubleVar(value=0.8 * Rs_var.get())

# Resolution
ttk.Label(input_frame, text="Resolution").grid(row=6, column=0, columnspan=2)
N_var = tk.IntVar(value=500)
ttk.Label(input_frame, text="N").grid(row=7, column=0)
ttk.Entry(input_frame, textvariable=N_var).grid(row=7, column=1)
N_l_var = tk.IntVar(value=12)
ttk.Label(input_frame, text="N_l").grid(row=8, column=0)
ttk.Entry(input_frame, textvariable=N_l_var).grid(row=8, column=1)
N_up_var = tk.IntVar(value=10)
ttk.Label(input_frame, text="N_up").grid(row=9, column=0)
ttk.Entry(input_frame, textvariable=N_up_var).grid(row=9, column=1)

# Trailing Edge Function
ttk.Label(input_frame, text="Trailing Edge Function").grid(row=10, column=0, columnspan=2)
trailing_edge_function_var = tk.StringVar(value="a + b*y**2 + c*y**4")
ttk.Entry(input_frame, textvariable=trailing_edge_function_var).grid(row=11, column=0, columnspan=2)

# Output file checkbutton
output_file_var = tk.BooleanVar()
ttk.Checkbutton(input_frame, text="Output file", variable=output_file_var).grid(row=12, column=0, columnspan=2)

# Run button
ttk.Button(input_frame, text="RUN", command=lambda: threading.Thread(target=run_script).start()).grid(row=13, column=0, columnspan=2)

# Create the image display frame
image_frame = ttk.Frame(root)
image_frame.grid(row=0, column=1, padx=10, pady=10, sticky="ne")
image_label = ttk.Label(image_frame)
image_label.grid(row=0, column=0)

# Create the output text frame
output_frame = ttk.Frame(root)
output_frame.grid(row=1, column=0, columnspan=2, padx=10, pady=10, sticky="s")
output_text = tk.Text(output_frame, height=10, width=80)
output_text.grid(row=0, column=0)

# Start the main loop
root.mainloop()
