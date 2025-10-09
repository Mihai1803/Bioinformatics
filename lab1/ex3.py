import tkinter as tk
from tkinter import filedialog, messagebox
from pathlib import Path
from collections import Counter

def read_fasta_sequence(path: str) -> str:
    seq_parts = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq_parts.append(line)
    return "".join(seq_parts).upper()

class FastaGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("FASTA Reader")
        self.geometry("640x480")

        # Top controls
        top = tk.Frame(self)
        top.pack(fill="x", padx=8, pady=8)

        self.path_var = tk.StringVar(value="No file selected")
        tk.Button(top, text="Open FASTA…", command=self.choose_file).pack(side="left")
        tk.Button(top, text="Analyze", command=self.analyze).pack(side="left", padx=6)
        tk.Button(top, text="Clear", command=self.clear_output).pack(side="left")

        tk.Label(top, textvariable=self.path_var, anchor="w").pack(side="left", padx=10)

        # Output box
        self.output = tk.Text(self, wrap="word", font=("Consolas", 10))
        self.output.pack(expand=True, fill="both", padx=8, pady=8)

        # Status bar
        self.status = tk.StringVar(value="Ready")
        tk.Label(self, textvariable=self.status, anchor="w", relief="sunken").pack(fill="x", side="bottom")

        self.file_path = None

    def choose_file(self):
        path = filedialog.askopenfilename(
            title="Select FASTA file",
            filetypes=[("FASTA files", "*.fasta *.fa *.txt"), ("All files", "*.*")]
        )
        if path:
            self.file_path = path
            self.path_var.set(path)
            self.status.set("File selected. Click Analyze.")

    def analyze(self):
        if not self.file_path or not Path(self.file_path).exists():
            messagebox.showwarning("No file", "Please select a valid FASTA file first.")
            return

        try:
            seq = read_fasta_sequence(self.file_path)
        except Exception as e:
            messagebox.showerror("Read error", f"Could not read file:\n{e}")
            return

        self.output.delete("1.0", "end")

        if not seq:
            self.output.insert("end", "No sequence found (empty file or only headers).\n")
            self.status.set("Done")
            return

        counts = Counter(seq)
        total = sum(counts.values())

        canonical = ["A", "C", "G", "T"]
        others = sorted([b for b in counts.keys() if b not in canonical])
        ordered = [b for b in canonical if b in counts] + others

        # Print results
        self.output.insert("end", "Alphabet (unique symbols):\n")
        self.output.insert("end", ", ".join(ordered) + "\n\n")

        self.output.insert("end", "Relative frequencies:\n")
        for b in ordered:
            frac = counts[b] / total if total else 0.0
            self.output.insert("end", f"{b}: {counts[b]} / {total} = {frac:.6f} | {frac*100:.3f}%\n")
        self.output.insert("end", "\n")

        self.status.set("Done")

    def clear_output(self):
        self.output.delete("1.0", "end")
        self.status.set("Cleared")

if __name__ == "__main__":
    FastaGUI().mainloop()
