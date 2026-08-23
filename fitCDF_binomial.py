import ROOT
import numpy as np
from scipy.optimize import minimize
from scipy.special import erf
import argparse
import csv

def gaussian_cdf(x, mu, sigma):

    return 0.5 * (
        1.0 + erf(
            (x - mu) /
            (np.sqrt(2.0) * sigma)
        )
    )


def binomial_nll(params, x, passed, total):

    mu, sigma = params

    if sigma <= 0:
        return 1e100

    p = gaussian_cdf(x, mu, sigma)

    p = np.clip(p, 1e-12, 1.0 - 1e-12)

    return -np.sum(
        passed * np.log(p)
        + (total - passed) * np.log(1.0 - p)
    )



def get_tefficiency_data(eff):

    h_total = eff.GetTotalHistogram()
    h_passed = eff.GetPassedHistogram()

    n_bins = h_total.GetNbinsX()

    x = []
    passed = []
    total = []
    efficiency = []
    err_low = []
    err_high = []

    for i in range(1, n_bins + 1):

        n = h_total.GetBinContent(i)

        if n <= 0:
            continue

        x.append(h_total.GetBinCenter(i))
        total.append(n)
        passed.append(h_passed.GetBinContent(i))

        efficiency.append(eff.GetEfficiency(i))
        err_low.append(eff.GetEfficiencyErrorLow(i))
        err_high.append(eff.GetEfficiencyErrorUp(i))
    return (
        np.asarray(x, dtype=float),
        np.asarray(passed, dtype=float),
        np.asarray(total, dtype=float),
    )

def fit_turnon(eff, mu_guess=None, sigma_guess=20.):

    x, passed, total = get_tefficiency_data(eff)

    if mu_guess is None:
        efficiency = passed / total

        idx = np.argmin(
            np.abs(efficiency - 0.5)
        )

        mu_guess = x[idx]

    result = minimize(
        binomial_nll,
        [mu_guess, sigma_guess],
        args=(x, passed, total),
        method="Nelder-Mead"
    )

    return result.x, result
    
def shift_tefficiency(eff, direction):

    h_total = eff.GetTotalHistogram()
    h_pass = eff.GetPassedHistogram()

    h_total_new = h_total.Clone()
    h_pass_new = h_pass.Clone()

    h_total_new.SetDirectory(0)
    h_pass_new.SetDirectory(0)

    n_bins = h_total.GetNbinsX()

    for i in range(1, n_bins + 1):

        total = h_total.GetBinContent(i)

        if total <= 0:
            continue

        efficiency = eff.GetEfficiency(i)

        if direction == "up":
            error = eff.GetEfficiencyErrorUp(i)
            shifted_efficiency = efficiency + error

        elif direction == "down":
            error = eff.GetEfficiencyErrorLow(i)
            shifted_efficiency = efficiency - error

        else:
            raise ValueError(
                "direction must be 'up' or 'down'"
            )

        shifted_efficiency = np.clip(
            shifted_efficiency,
            0.0,
            1.0
        )

        passed_new = int(
            np.rint(
                shifted_efficiency * total
            )
        )

        h_pass_new.SetBinContent(
            i,
            passed_new
        )

    return ROOT.TEfficiency(
        h_pass_new,
        h_total_new
    )
    
def plot_turnon(fits):

    canvas = ROOT.TCanvas(
        "c_turnon",
        "Trigger Turn-On",
        800,
        600
    )

    canvas.SetLeftMargin(0.12)
    canvas.SetBottomMargin(0.12)

    # --------------------------------------------------------
    # Get x range from first TEfficiency
    # --------------------------------------------------------

    eff0 = fits[0][0]

    h_total = eff0.GetTotalHistogram()

    xmin = h_total.GetXaxis().GetXmin()
    xmax = h_total.GetXaxis().GetXmax()

    # --------------------------------------------------------
    # Frame
    # --------------------------------------------------------

    frame = canvas.DrawFrame(
        xmin,
        0.0,
        xmax,
        1.05
    )

    frame.SetTitle("")
    frame.GetXaxis().SetTitle("MET [GeV]")
    frame.GetYaxis().SetTitle("Trigger efficiency")

    # --------------------------------------------------------
    # Keep ROOT/Python objects alive
    # --------------------------------------------------------

    canvas._turnon_objects = [frame]

    # --------------------------------------------------------
    # Legend
    # --------------------------------------------------------

    legend = ROOT.TLegend(
        0.50,
        0.18,
        0.88,
        0.42
    )

    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    # --------------------------------------------------------
    # Draw each TEfficiency + fitted CDF
    # --------------------------------------------------------

    for eff, fit, label, color, linestyle in fits:

        # --------------------------------------------
        # Draw efficiency
        # --------------------------------------------

        eff.SetMarkerStyle(20)
        eff.SetMarkerSize(0.8)

        eff.SetMarkerColor(color)
        eff.SetLineColor(color)

        eff.Draw("P SAME")

        # --------------------------------------------
        # Fit parameters
        # --------------------------------------------

        mu, sigma = fit

        # --------------------------------------------
        # Make smooth CDF
        # --------------------------------------------

        x = np.linspace(
            xmin,
            xmax,
            500
        )

        y = gaussian_cdf(
            x,
            float(mu),
            float(sigma)
        )

        graph = ROOT.TGraph(
            len(x),
            np.asarray(x, dtype="double"),
            np.asarray(y, dtype="double")
        )

        graph.SetLineColor(color)
        graph.SetLineWidth(3)
        graph.SetLineStyle(linestyle)

        graph.Draw("L SAME")

        # --------------------------------------------
        # Keep objects alive!
        # --------------------------------------------

        canvas._turnon_objects.append(eff)
        canvas._turnon_objects.append(graph)

        # --------------------------------------------
        # Legend
        # --------------------------------------------

        legend.AddEntry(
            graph,
            label,
            "l"
        )

    legend.Draw()

    canvas._turnon_objects.append(legend)

    canvas.Modified()
    canvas.Update()

    return canvas
    
def ProcessEfficiency( eff, writer ):
	
	# ------------------------------------------------------------
	# Nominal fit
	# ------------------------------------------------------------

	fit_nominal, _ = fit_turnon(eff)


	# ------------------------------------------------------------
	# Up/down variations
	# ------------------------------------------------------------

	eff_up = shift_tefficiency(
    	eff,
    	"up"
	)

	eff_down = shift_tefficiency(
    	eff,
    	"down"	
	)


	# ------------------------------------------------------------
	# Fit variations
	# ------------------------------------------------------------

	fit_up, _ = fit_turnon(eff_up)

	fit_down, _ = fit_turnon(eff_down)


	# ------------------------------------------------------------
	# Plot everything
	# ------------------------------------------------------------

	fits = [
  	  (
     	   eff,
    	    fit_nominal,
        	"Nominal",
        	ROOT.kBlack,
        	1
    	),

    	(
        	eff_up,
        	fit_up,
        	"Up",
        	ROOT.kRed,
        	2
    	),

    	(
        	eff_down,
        	fit_down,
        	"Down",
        	ROOT.kBlue,
        	2
   		 )
	]

	canvas = plot_turnon(fits)
	canvas.SaveAs(eff.GetName()+".pdf")
	csveffname = eff.GetName()
	HTlow = "-1"
	HThigh = "-1"
	year = "-1"
	mu_nom,sigma_nom = fit_nominal
	mu_up,sigma_up = fit_up
	mu_down,sigma_down = fit_down
	rownom = [csveffname, HTlow, HThigh, year,str(mu_nom), str(sigma_nom), "nom"] 
	rowup = [csveffname, HTlow, HThigh, year, str(mu_up), str(sigma_up), "up"]
	rowdown = [csveffname, HTlow, HThigh, year, str(mu_down), str(sigma_down), "down"]
	
	writer.writerow(rownom)
	writer.writerow(rowup)
	writer.writerow(rowdown)


parser = argparse.ArgumentParser(description="Description of what your script does.")

    # 2. Add a required positional argument
parser.add_argument("-i", "--input", type=str, help="Path to the input ROOT file CSV (ROOTFile,TEfficiencyName)")

    # 3. Add an optional value argument (with a default value)
parser.add_argument("-o", "--output", type=str, default="output.txt", help="Path to the output csv name") 

    # 6. Parse the arguments
args = parser.parse_args()

    # 7. Use the arguments in your script
print(f"Processing: {args.input}")
print(f"Output destination: {args.output}")


with open(args.input, 'r', newline='') as infile, open(args.output, 'w', newline='') as outfile:
    reader = csv.reader(infile)
    writer = csv.writer(outfile)
    
    # Optional: Read and write the header row
    header=["triggerName","HTlow","HThigh","year","mu","sigma","syst"]
    writer.writerow(header)
    
    # Create a reader object
    reader = csv.reader(infile)
    
    # Optional: Skip the header row if you don't want to print it
    header = next(reader)
    print(f"Headers: {header}")
    
    # Loop through each row
    for row in reader:
        print(row)
        tfile = ROOT.TFile.Open(row[0])
        eff = tfile.Get(row[1])
        ProcessEfficiency( eff, writer ) 

  


