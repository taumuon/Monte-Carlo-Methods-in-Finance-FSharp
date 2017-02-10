// This code is based on matlab code provided through the course "Monte Carlo Methods in Finance".
// https://iversity.org/my/courses/monte-carlo-methods-in-finance/
// and Olaf Smits's Python conversion 
// http://nbviewer.ipython.org/github/olafSmits/MonteCarloMethodsInFinance/blob/master/Week%201.ipynb?create=1

open System
open Deedle
open FSharp.Charting

let readFrame (stock:string) =
    let frame = Frame.ReadCsv(stock + ".csv") |> Frame.indexRowsDate "Date" |> Frame.orderRows
    let newFrame = frame.Columns.[ ["Date"; "Adj Close" ] ]
    let logRet = newFrame.GetSeries<float>("Adj Close")
                 |> Series.pairwiseWith(fun k (v1, v2) -> Math.Log(v2 / v1))
    newFrame.AddSeries("logRet", logRet)
    newFrame.RenameSeries(fun s -> (stock + s).Replace(" ", ""))
    newFrame

[<EntryPoint>]
let main argv = 
    let frame = readFrame("IBM").Join(readFrame("GOOG"), JoinKind.Inner)
                                .Join(readFrame("SI"), JoinKind.Inner)
    frame.DropSeries("GOOGDate")
    frame.DropSeries("SIDate")

    let mean = frame?IBMAdjClose |> Series.mean

    let chart = frame.GetAllSeries() 
                |> Seq.map( fun (KeyValue(k,v)) -> Chart.Line(v |> Series.observations, Name=k))
                |> Chart.Combine
                |> Chart.WithLegend (Docking = ChartTypes.Docking.Top)
    chart.ShowChart()

    let logRetChart = frame.GetAllSeries()
                            |> Seq.filter (fun (KeyValue(k,v)) -> k.EndsWith("logRet"))
                            |> Seq.map( fun (KeyValue(k,v)) -> Chart.Line(v |> Series.observations, Name=k))
                            |> Chart.Combine
                            |> Chart.WithLegend (Docking = ChartTypes.Docking.Top)
    logRetChart.ShowChart()

    printfn "%A" mean
    System.Windows.Forms.Application.Run()
    0