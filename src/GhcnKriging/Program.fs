open MathNet.Numerics.LinearAlgebra

open System.Runtime.InteropServices
open System.Collections.Immutable
open System.IO

open Angara
open Angara.Charting
open Angara.Data
open Angara.Statistics
open Angara.Filzbach
open Angara.Serialization

 type KrigingResult = {Value: float[]; Error: float[]}

[<EntryPoint>]
let main argv = 
    // Preparing data
    
    // Reading a CSV data file
    let ghcn = Table.Load("GHCN_monthly.csv")
    // Filtering data and Fahrenheit/Celsium conversion
    let latMin = 20.0
    let latMax = 50.0
    let lonMin = 0.0
    let lonMax = 30.0
    let ghcn_filtered = ghcn |> Table.Filter ["lat"; "lon"; "temp_mean"] (fun lat lon temp -> temp <> -9999.0 && latMin <= lat && lat <= latMax && lonMin <= lon && lon <= lonMax) 
                        |> Table.MapToColumn "temp_mean" [|"temp_mean"|] (fun f -> (f - 32.0)/1.8)
    let lat = ghcn_filtered.["lat"].Rows.AsReal |> Seq.toArray
    let lon = ghcn_filtered.["lon"].Rows.AsReal|> Seq.toArray
    let temp = ghcn_filtered.["temp_mean"].Rows.AsReal |>Seq.toArray
    // Forming lat and lon grid for kriging
    let lat_count = 50
    let lon_count = 50
    let latGrid1 = [| for i in 0..(lat_count - 1)  -> latMin + (latMax - latMin) * double(i) / double(lat_count - 1) |]
    let lonGrid1 = [| for i in 0..(lon_count - 1)-> lonMin + (lonMax - lonMin) * double(i) / double(lon_count - 1) |]
    let latGrid = [| for i in 0..(lat_count * lon_count - 1) -> latMin + (latMax - latMin) * double(i % lat_count) / double(lat_count - 1) |]
    let lonGrid = [| for i in 0..(lat_count * lon_count - 1) -> lonMin + (lonMax - lonMin) * double(i / lon_count) / double(lon_count - 1) |]
    // Visualizaton of processed data as markers chart
    let markers_plots =  [ 
                                Plot.markers(lon, lat, 
                                    color = MarkersColor.Values temp, colorPalette = "Blue,Green,Red", 
                                    shape = MarkersShape.Circle, displayName = "lon/lat/temp",
                                    titles = Titles.markers("lon", "lat"))
                         ]
    let markers_chart = Chart.ofList markers_plots
    Angara.Base.Init();
    Angara.Html.Save "markers chart.html" markers_chart
    System.Diagnostics.Process.Start("markers chart.html") |> ignore

    //Building SVGram
    let maxDist = 2000.0
    let binWidth = 100.0

    let sphdist(lat1 : double, lon1 : double, lat2 : double, lon2 : double) = 
        let deg2rad a = 3.1415 * a / 180.0
        let lat1 = deg2rad lat1
        let lat2 = deg2rad lat2
        let lon1 = deg2rad lon1
        let lon2 = deg2rad lon2
        let dlat = abs(lat1 - lat2)
        let dlon = abs(lon1 - lon2)
        let slat = sin(dlat/2.0)
        let slon = sin(dlon/2.0)
        let dsigma = 2.0 * asin(sqrt(slat * slat + cos(lat1) * cos(lat2) * slon * slon))
        dsigma * 6371.0

    let expGamma n s r h = 
        let e = (s - n) * (1.0 - exp (-3.0*h/r))
        if h > 0.0 then e + n else e

    let buildSVGram (lat  : double[]) (lon : double[]) (temp : double[]) = 
        let n = temp.Length
        let count = Array.zeroCreate<int> (int(maxDist / binWidth))
        let sum = Array.zeroCreate<double> (int (maxDist / binWidth))
        for i in 0..n - 1 do
            for j in i + 1 .. n - 1 do
                let dist = int(round(sphdist(lat.[i],lon.[i],lat.[j],lon.[j]) / binWidth))
                if dist < count.Length then
                    sum.[dist] <- sum.[dist] + sqrt((temp.[i] - temp.[j])*(temp.[i] - temp.[j]))/2.0
                    count.[dist] <- count.[dist] + 1
        Array.zip count sum 
        |> Array.filter(fun (c,s) -> c > 10)
        |> Array.mapi(fun i (c,s) -> double(binWidth) / 2.0 + double(i) * binWidth, s / double(c))
        |> Array.unzip

    let dist, gamma = buildSVGram lat lon temp

    let expLogLikehood (values:Parameters) =
        Array.zip dist gamma 
        |> Array.map(fun (dist,gamma) -> 
                         let gamma2 = expGamma (values.GetValue "n") (values.GetValue "s") (values.GetValue "r") dist
                         let stdev = values.GetValue "sigma"
                         if System.Double.IsNaN(gamma2) 
                         then 0.0 
                         else let dev = (gamma-gamma2)/stdev in -0.5*dev*dev - log(sqrt(2.0*System.Math.PI)*stdev))
        |> Array.sum
    
    let RunFilzbach () =
            //Entering parameters for Filzbach
            let p0 = Parameters.Empty
            let p1 = p0.Add("n", Distribution.Uniform(0.0, 1000.0)) 
            let p2 = p1.Add("s", Distribution.Uniform(0.0, 1000.0))
            let p3 = p2.Add("r",  Distribution.Uniform(0.0, 1000.0)) 
            let p4 = p3.Add("sigma", Distribution.Uniform(0.0, 100.0)) 
            //Running Filzbach 
            let run = Sampler.runmcmc(p4, expLogLikehood, 1000, 5000) 

            //Retrieving Filzbach results
            let getsummary (result : SamplerResult) (parameter_name : string) =
                let parameter_index =  (run.sampler.Parameters.GetDefinition parameter_name).index
                qsummary (result.samples |> Seq.map (fun {values=sample} -> sample.[parameter_index]))

            ((getsummary run "n").median, (getsummary run "s").median, (getsummary run "r").median) 

    // Running Filzbach or restoring results from checkpoint with Reinstate
    let checkpoint_name = "Filzbach parameters"
    let n, s ,r = Angara.ReinstateServices.Reinstate checkpoint_name RunFilzbach
    // Generate exponential model 
    let dist2, gamma2 = [| for i in 0..int(maxDist) -> float(i), expGamma n s r (double(i)) |] |> Array.unzip
    // Visualization of svgram data and exponential model
    let semivariogram_plots =  [ 
                                 Plot.markers(MarkersX.Values dist, MarkersY.Values gamma, displayName = "gamma") ;
                                 Plot.line(dist2, gamma2, displayName = "gamma2") 
                               ]
    let semivariogram_chart = Chart.ofList semivariogram_plots
    Angara.Html.Save "semivariogram chart.html" semivariogram_chart
    System.Diagnostics.Process.Start("semivariogram chart.html") |> ignore

    //Kriging
    let krige (xi : double[], yi : double[], x : double[], y : double[], f : double[], gamma: double -> double, dist: double * double * double * double -> double) =
        let N = x.Length
        let Г1_column i = Array.init (N+1) (fun j -> if i >= N then if j >= N then 0.0 else 1.0
                                                               else if j >= N then 1.0 else -gamma(dist(x.[j], y.[j], x.[i], y.[i])))

        let Г1_inv = CreateMatrix.DenseOfColumnArrays(Array.init (N+1) Г1_column).Inverse()
        let Г0 i j = CreateVector.DenseOfArray(Array.init (N+1) (fun k -> if k >= N then 1.0 else -gamma(dist(xi.[i], yi.[j], x.[k], y.[k]))))
        let K = xi.Length
        let run = Array2D.init (xi.Length) (yi.Length) (fun i j -> let v = Г0 i j
                                                                   let w = Г1_inv.Multiply(v)
                                                                   let mutable f0 = 0.0
                                                                   let mutable f1 = 0.0
                                                                   for k in 0..N-1 do
                                                                                f0 <- f0 + w.At(k) * f.[k]
                                                                                f1 <- f1 + w.At(k) * (v).[k]
                                                                   (f0, -f1)) |> Seq.cast<float*float> 
        let res = run|> Seq.map (fun (a,b) -> a) |> Seq.toArray
        let err = run|> Seq.map (fun (a,b) -> b) |> Seq.toArray
        {Value = res; Error = err}

    // Running Simple Kriging
    let Kriege = krige(lonGrid1, latGrid1, lon, lat, temp, expGamma n s r, sphdist)
    let grid = Kriege.Value
    let error =  Kriege.Error

    // Forming quantiles
    let heatmap_quantiles = 
        { median = grid
          lower68 = Array.init (grid.Length) (fun i -> grid.[i] - error.[i])
          upper68 = Array.init (grid.Length) (fun i -> grid.[i] + error.[i])
          lower95 = Array.init (grid.Length) (fun i -> grid.[i] - 3.0 * error.[i])
          upper95 = Array.init (grid.Length) (fun i -> grid.[i] + 3.0 * error.[i]) }
    // Visualizing kriging result and variation
    let heatmap_chart = [
                            [  
                                Plot.markers(lon, lat, color = MarkersColor.Values temp, colorPalette = "Blue,Green,Red", shape = MarkersShape.Circle, displayName = "lon/lat/temp", 
                                             titles = Titles.markers("lon","lat","temp"));
                                Plot.heatmap(lonGrid, latGrid, HeatmapValues.TabularUncertainValues heatmap_quantiles, colorPalette = "Blue,Green,Red", titles = Titles.heatmap("lon", "lat"))
                            ] |> Chart.ofList;
                            [ 
                                Plot.markers(MarkersX.Values lon, MarkersY.Values lat, shape = MarkersShape.Cross, titles = Titles.markers("lon", "lat", "temp"));
                                Plot.heatmap(lonGrid, latGrid, error, colorPalette = "Blue,Green,Red", titles = Titles.heatmap("lon", "lat"))
                            ] |> Chart.ofList
                        ]            
    Angara.Html.Save "heatmap chart.html" heatmap_chart
    System.Diagnostics.Process.Start("heatmap chart.html") |> ignore

    System.Console.ReadLine() |> ignore
    0 // return an integer exit code
