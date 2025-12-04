import tkinter as tk
from tkinter import ttk, filedialog, messagebox

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
import matplotlib as mpl

mpl.rcParams.update({
    "font.size": 8,
})


def clean_seq(seq: str) -> str:
    return ''.join(b for b in seq.upper() if b in "ACGT")


def cg_total(seq: str) -> float:
    s = clean_seq(seq)
    if not s:
        return 0.0
    cg = sum(1 for b in s if b in "CG")
    return 100.0 * cg / len(s)


def cpg_total(seq: str) -> int:
    s = clean_seq(seq)
    return sum(1 for i in range(len(s) - 1) if s[i:i+2] == "CG")


def cg_window_relative(window: str, total_len: int) -> float:
    w = clean_seq(window)
    if not w or total_len == 0:
        return 0.0
    cg = sum(1 for b in w if b in "CG")
    return 100.0 * cg / total_len


def cpg_window(window: str) -> int:
    w = clean_seq(window)
    return sum(1 for i in range(len(w) - 1) if w[i:i+2] == "CG")


def kappa_ic_window(window: str) -> float:
    A = clean_seq(window)
    L = len(A)
    if L < 2:
        return 0.0

    N = L - 1
    T = 0.0
    for u in range(1, N + 1):
        lengthB = L - u
        C = 0
        for i in range(lengthB):
            if A[i] == A[i + u]:
                C += 1
        T += (C / lengthB) * 100.0
    return T / N


def sliding_metrics(seq: str, window_len: int, step: int = 1):
    s = clean_seq(seq)
    if len(s) < window_len:
        raise ValueError("Sequence shorter than window length.")

    total_len = len(s)

    positions, cg_list, ic_list, cpg_list = [], [], [], []

    for start in range(0, len(s) - window_len + 1, step):
        w = s[start:start + window_len]
        center_pos = start + window_len // 2 + 1

        positions.append(center_pos)
        cg_list.append(cg_window_relative(w, total_len))
        ic_list.append(kappa_ic_window(w))
        cpg_list.append(cpg_window(w))

    return positions, cg_list, ic_list, cpg_list


def center_of_weight(xs, ys):
    if not xs:
        return (0.0, 0.0)
    x_c = sum(xs) / len(xs)
    y_c = sum(ys) / len(ys)
    return (x_c, y_c)


class PromKappaGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Mini PromKappa - DNA promoter patterns")
        self.geometry("1200x700")

        self.centers_cg = []
        self.centers_cpg = []

        self.create_widgets()
        self.create_plots()

    def create_widgets(self):
        left = ttk.Frame(self, padding=5)
        left.pack(side=tk.LEFT, fill=tk.Y)

        params = ttk.LabelFrame(left, text="Parameters", padding=5)
        params.pack(side=tk.TOP, fill=tk.X, pady=5)

        ttk.Label(params, text="Window length:").grid(row=0, column=0, sticky="w")
        self.window_len_var = tk.StringVar(value="30")
        ttk.Entry(params, textvariable=self.window_len_var, width=6).grid(
            row=0, column=1, sticky="w"
        )

        ttk.Label(params, text="Window step:").grid(row=1, column=0, sticky="w")
        self.window_step_var = tk.StringVar(value="1")
        ttk.Entry(params, textvariable=self.window_step_var, width=6).grid(
            row=1, column=1, sticky="w"
        )

        seq_frame = ttk.LabelFrame(left, text="Promoter sequence", padding=5)
        seq_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True, pady=5)

        ttk.Label(seq_frame, text="DNA sequence (A/C/G/T):").pack(anchor="w")
        self.seq_text = tk.Text(seq_frame, width=40, height=18, wrap="word")
        self.seq_text.pack(fill=tk.BOTH, expand=True)

        example_S = (
            "CGGACTGATCTATCTAAAAAAAAAAAAAAAAAAAAAAAAAAACGTAGCATCTATCGAT"
            "CTATCTAGCGATCTATCTACTACG"
        )
        self.seq_text.insert("1.0", example_S)

        btn_frame = ttk.Frame(seq_frame)
        btn_frame.pack(fill=tk.X, pady=3)

        ttk.Button(btn_frame, text="Load FASTA...", command=self.load_fasta).pack(
            side=tk.LEFT, padx=2
        )
        ttk.Button(btn_frame, text="Analyze", command=self.analyze_current).pack(
            side=tk.LEFT, padx=2
        )

        res_frame = ttk.LabelFrame(left, text="Results (current promoter)", padding=5)
        res_frame.pack(side=tk.TOP, fill=tk.X, pady=5)

        self.cg_label = ttk.Label(res_frame, text="(C+G)% (whole) : ---")
        self.cg_label.pack(anchor="w")
        self.cpg_label = ttk.Label(res_frame, text="CpG (CG) count (whole): ---")
        self.cpg_label.pack(anchor="w")
        self.ic_label = ttk.Label(res_frame, text="Kappa IC (whole): ---")
        self.ic_label.pack(anchor="w")
        self.center_cg_label = ttk.Label(res_frame, text="Center (C+G%, IC): ---")
        self.center_cg_label.pack(anchor="w")
        self.center_cpg_label = ttk.Label(res_frame, text="Center (CG, IC): ---")
        self.center_cpg_label.pack(anchor="w")

        ttk.Button(
            res_frame, text="Clear centers chart", command=self.clear_centers
        ).pack(pady=(5, 0), anchor="w")

    def create_plots(self):
        mid = ttk.Frame(self, padding=5)
        mid.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self.fig = Figure(figsize=(9, 6), dpi=100)

        gs = GridSpec(3, 2, height_ratios=[1.2, 1, 1], figure=self.fig)

        self.fig.subplots_adjust(
            left=0.08,
            right=0.99,
            top=0.96,
            bottom=0.07,
            wspace=0.40,
            hspace=0.55,
        )

        self.ax_top = self.fig.add_subplot(gs[0, :])
        self.ax_p11 = self.fig.add_subplot(gs[1, 0])
        self.ax_p12 = self.fig.add_subplot(gs[1, 1])
        self.ax_p21 = self.fig.add_subplot(gs[2, 0])
        self.ax_p22 = self.fig.add_subplot(gs[2, 1])

        for ax in (self.ax_top, self.ax_p11, self.ax_p12, self.ax_p21, self.ax_p22):
            ax.tick_params(labelsize=7)

        self.canvas = FigureCanvasTkAgg(self.fig, master=mid)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.canvas.draw()

    def load_fasta(self):
        path = filedialog.askopenfilename(
            title="Open FASTA file",
            filetypes=[("FASTA files", "*.fa *.fasta *.txt"), ("All files", "*.*")],
        )
        if not path:
            return

        try:
            with open(path, "r") as f:
                lines = f.readlines()
        except Exception as e:
            messagebox.showerror("Error", f"Could not read file:\n{e}")
            return

        seq_parts = []
        for line in lines:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq_parts.append(line)

        if not seq_parts:
            messagebox.showwarning("Warning", "No sequence found in FASTA.")
            return

        seq = "".join(seq_parts)
        self.seq_text.delete("1.0", tk.END)
        self.seq_text.insert(tk.END, seq)

    def clear_centers(self):
        self.centers_cg.clear()
        self.centers_cpg.clear()
        self.ax_p21.clear()
        self.ax_p22.clear()
        self.ax_p21.set_title("Centers (C+G%, IC)", fontsize=9)
        self.ax_p22.set_title("Centers (CG, IC)", fontsize=9)
        self.ax_p21.set_xlabel("C+G%", fontsize=8)
        self.ax_p21.set_ylabel("Kappa IC", fontsize=8)
        self.ax_p22.set_xlabel("CG", fontsize=8)
        self.ax_p22.set_ylabel("Kappa IC", fontsize=8)
        self.ax_p21.grid(True)
        self.ax_p22.grid(True)
        self.ax_p21.tick_params(labelsize=7)
        self.ax_p22.tick_params(labelsize=7)
        self.canvas.draw()

    def analyze_current(self):
        seq = self.seq_text.get("1.0", tk.END).strip()
        if not seq:
            messagebox.showwarning("Warning", "Please enter a DNA sequence.")
            return

        try:
            win_len = int(self.window_len_var.get())
            step = int(self.window_step_var.get())
            if win_len <= 1 or step <= 0:
                raise ValueError
        except ValueError:
            messagebox.showerror(
                "Error", "Window length must be >1 and step must be >0 (integers)."
            )
            return

        cg_tot = cg_total(seq)
        cpg_tot = cpg_total(seq)
        ic_tot = kappa_ic_window(seq)

        self.cg_label.config(text=f"(C+G)% (whole) : {cg_tot:.2f}")
        self.cpg_label.config(text=f"CpG (CG) count (whole): {cpg_tot}")
        self.ic_label.config(text=f"Kappa IC (whole): {ic_tot:.2f}")

        try:
            positions, cg_list, ic_list, cpg_list = sliding_metrics(seq, win_len, step)
        except ValueError as e:
            messagebox.showerror("Error", str(e))
            return

        cx_cg, cy_cg = center_of_weight(cg_list, ic_list)
        cx_cpg, cy_cpg = center_of_weight(cpg_list, ic_list)

        self.center_cg_label.config(
            text=f"Center (C+G%, IC): ({cx_cg:.2f}, {cy_cg:.2f})"
        )
        self.center_cpg_label.config(
            text=f"Center (CG, IC): ({cx_cpg:.2f}, {cy_cpg:.2f})"
        )

        self.centers_cg.append((cx_cg, cy_cg))
        self.centers_cpg.append((cx_cpg, cy_cpg))

        self.ax_top.clear()
        self.ax_top.plot(positions, ic_list, label="Kappa IC")
        self.ax_top.plot(positions, cg_list, label="(C+G)% (relative)")
        self.ax_top.plot(positions, cpg_list, label="CpG (CG) count")
        self.ax_top.set_title("Sliding-window metrics", fontsize=10)
        self.ax_top.set_xlabel("Position (bp)", fontsize=8)
        self.ax_top.set_ylabel("Value", fontsize=8)
        self.ax_top.grid(True)
        self.ax_top.legend(loc="upper right", fontsize=7)
        self.ax_top.tick_params(labelsize=7)

        def draw_pattern(ax, xs, ys, xlabel, center_x, center_y, color, title):
            ax.clear()
            ax.scatter(xs, ys, s=6, color=color)
            ax.axvline(center_x, linestyle="--", linewidth=0.7)
            ax.axhline(center_y, linestyle="--", linewidth=0.7)
            ax.set_title(title, fontsize=9)
            ax.set_xlabel(xlabel, fontsize=8)
            ax.set_ylabel("Kappa IC", fontsize=8)
            ax.grid(True)
            ax.tick_params(labelsize=7)

        draw_pattern(
            self.ax_p11, cg_list, ic_list,
            "C+G% (relative)", cx_cg, cy_cg, "tab:blue",
            "Pattern: IC vs C+G%"
        )

        draw_pattern(
            self.ax_p12, cpg_list, ic_list,
            "CG (count)", cx_cpg, cy_cpg, "tab:red",
            "Pattern: IC vs CG count"
        )

        def draw_centers(ax, centers, xlabel, title):
            ax.clear()
            if centers:
                xs = [c[0] for c in centers]
                ys = [c[1] for c in centers]
                ax.scatter(xs, ys, s=25, marker="o")
            ax.set_title(title, fontsize=9)
            ax.set_xlabel(xlabel, fontsize=8)
            ax.set_ylabel("Kappa IC", fontsize=8)
            ax.grid(True)
            ax.tick_params(labelsize=7)

        draw_centers(
            self.ax_p21, self.centers_cg,
            "C+G% (center of each pattern)", "Centers (C+G%, IC)"
        )

        draw_centers(
            self.ax_p22, self.centers_cpg,
            "CG (center of each pattern)", "Centers (CG, IC)"
        )

        self.canvas.draw()


if __name__ == "__main__":
    app = PromKappaGUI()
    app.mainloop()
