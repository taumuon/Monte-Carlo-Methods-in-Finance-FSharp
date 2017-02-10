open MathNet.Numerics.Distributions
open MathNet.Numerics.Statistics

let callPayoff strike price = max (price - strike) 0.0

let europeanPayoff payoff assetPath = assetPath |> Seq.last |> payoff

let europeanCallPayoff strike assetPath = assetPath |> europeanPayoff (callPayoff strike)

let asianArithmeticMeanCallPayoff strike (assetPath:seq<float>) = assetPath.Mean() |> callPayoff strike

let upAndOutPayoff payoff upperBarrier assetPath =
    if assetPath |> Seq.exists (fun x -> x > upperBarrier) then 0.0 else payoff assetPath

let upAndOutEuropeanCallPayoff strike = upAndOutPayoff (europeanCallPayoff strike)

let doubleBarrierPayoff payoff upperBarrier lowerBarrier assetPath =
    if assetPath |> Seq.exists (fun x -> x > upperBarrier || x < lowerBarrier) then 0.0 else payoff assetPath

let doubleBarrierEuropeanCallPayoff strike = doubleBarrierPayoff (europeanCallPayoff strike)

let getAssetPath S0 r deltaT sigma (normal:Normal) numSamples =
    Seq.unfold (fun S -> let sNew = (S * exp(((r - (0.5 * sigma * sigma)) * deltaT) + (sigma * sqrt(deltaT) * normal.Sample())))
                         Some(sNew, sNew)) S0
    |> Seq.take numSamples

let phi x = 
    let normal = new Normal()
    normal.CumulativeDistribution(x)

// Exact option pricing:

let priceEuropeanCall S0 strike r T sigma =
    let discountedK = strike * exp(-r * T)
    let totalVolatility = sigma * sqrt(T)
    let d_minus = log(S0 / discountedK) / totalVolatility
    let d_plus = d_minus + (0.5 * totalVolatility)
    let d_minus2 = d_minus - (0.5 * totalVolatility)
    (S0 * phi(d_plus)) - (discountedK * phi(d_minus2))

let computeD S0 strike r sigma T =
    let d1 = (log (S0 / strike) + ((r + (sigma * sigma / 2.0)) * T)) / (sigma * sqrt(T))
    let d0 = d1 - (sigma * sqrt(T))
    (d0, d1)

let priceBarrierUpAndOutCall S0 strike B r T sigma =
    let normal = new Normal()
    let mu = r - (0.5 * sigma * sigma)
    let nu = 2.0 * (mu / (sigma * sigma))
    let aux1 = (B / S0) ** (nu + 2.0)
    let aux2 = (B / S0) ** nu
    let d = computeD S0 strike r sigma T
    let price = S0 * normal.CumulativeDistribution(snd d) - (strike * exp(-r * T) * normal.CumulativeDistribution(fst d))
    let d2 = computeD S0 B r sigma T
    let price2 = price - (S0 * normal.CumulativeDistribution(snd d2)) + (strike * exp(-r * T) * normal.CumulativeDistribution(fst d2))
    let d3 = computeD (1.0/S0) (1.0/B) r sigma T
    let price3 = price2 + (aux1 * S0 * normal.CumulativeDistribution(snd d3)) - (aux2 * strike * exp(-r * T) * normal.CumulativeDistribution(fst d3))
    let d4 = computeD (1.0 / (strike * S0))  (1.0 / (B * B)) r sigma T
    let price4 = price3 - (aux1 * S0 * normal.CumulativeDistribution(snd d4)) + (aux2 * strike * exp(-r * T) * normal.CumulativeDistribution(fst d4))
    price4

// Monte Carlo pricing

let priceAsianArithmeticMeanMC S0 strike r T sigma numTrajectories numSamples =
    let normal = new Normal(0.0, 1.0)
    let deltaT = T / float numSamples
    let payoffs = seq { for n in 1 .. numTrajectories do 
                        let assetPath = getAssetPath S0 r deltaT sigma normal numSamples |> Seq.toList
                        yield assetPath |> asianArithmeticMeanCallPayoff strike
                      }
    let discountFactor = exp(-r * T)
    let priceMC = discountFactor * payoffs.Mean()
    let stddevMC = discountFactor * payoffs.StandardDeviation() / sqrt(float numTrajectories)
    (priceMC, stddevMC)

let priceDoubleBarrierCallMC S0 strike upperBarrier lowerBarrier r T sigma numTrajectories numSamples =
    let normal = new Normal(0.0, 1.0)
    let deltaT = T / float numSamples
    let payoffs = seq { for n in 1 .. numTrajectories do 
                        let assetPath = getAssetPath S0 r deltaT sigma normal numSamples |> Seq.toList
                        let price = assetPath |> doubleBarrierEuropeanCallPayoff strike upperBarrier lowerBarrier
                        yield price
                      }
    let discountFactor = exp(-r * T)
    let priceMC = discountFactor * payoffs.Mean()
    let stddevMC = discountFactor * payoffs.StandardDeviation() / sqrt(float numTrajectories)
    (priceMC, stddevMC)

let controlVariate numTrajectories discountFactor (samplesFactory : unit -> seq<float>) controlExact (payoff : seq<float> -> float) (controlPayoff) =
    let samplePayoffs = seq { for n in 1 .. numTrajectories do 
                              let assetPath = samplesFactory() |> Seq.toList
                              let payoff = (payoff assetPath, controlPayoff assetPath)
                              yield payoff
                            }
                        |> Seq.toList
    let payoffs = samplePayoffs |> Seq.map (fun i -> fst i) |> Seq.toArray
    let payoffsControlVariate = samplePayoffs |> Seq.map (fun i -> snd i) |> Seq.toArray

    let variance = ArrayStatistics.PopulationCovariance(payoffs, payoffs)
    let varianceControlVariate = ArrayStatistics.PopulationCovariance(payoffsControlVariate, payoffsControlVariate)
    let covariance = ArrayStatistics.PopulationCovariance(payoffs, payoffsControlVariate)
    let correlation = covariance / sqrt(varianceControlVariate * variance)

    let priceMC = discountFactor * payoffs.Mean()
    let stddevMC = discountFactor * payoffs.StandardDeviation() / sqrt(float numTrajectories)

    let priceControlVariateMC = discountFactor * payoffsControlVariate.Mean()
    let stddevControlVariateMC = discountFactor * payoffsControlVariate.StandardDeviation() / sqrt(float numTrajectories)

    let price = priceMC - ((covariance / varianceControlVariate) * (priceControlVariateMC - controlExact))
    let stddev = stddevMC * sqrt(1.0 - (correlation * correlation))

    (price, stddev)

let priceDoubleBarrierCall_MC_CV_EU S0 strike upperBarrier lowerBarrier r T sigma numTrajectories numSamples =
    let normal = new Normal(0.0, 1.0)
    let discountFactor = exp(-r * T)
    let deltaT = T / float numSamples;

    let priceEuropeanExact = priceEuropeanCall S0 strike r T sigma
    let europeanCallfn = (fun assetPath -> europeanCallPayoff strike assetPath)

    let barrierCallPayoff = (fun assetPath -> doubleBarrierEuropeanCallPayoff strike upperBarrier lowerBarrier assetPath)
    controlVariate numTrajectories
                   discountFactor
                   (fun () -> getAssetPath S0 r deltaT sigma normal numSamples)
                   priceEuropeanExact
                   barrierCallPayoff
                   europeanCallfn

let priceDoubleBarrierCall_MC_CV_UpAndOut S0 strike upperBarrier lowerBarrier r T sigma numTrajectories numSamples =
    let normal = new Normal(0.0, 1.0)
    let discountFactor = exp(-r * T)
    let deltaT = T / float numSamples;

    let priceUpAndOutExact = priceBarrierUpAndOutCall S0 strike upperBarrier r T sigma

    let barrierCallPayoff = (fun assetPath -> doubleBarrierEuropeanCallPayoff strike upperBarrier lowerBarrier assetPath)
    let upAndOut = (fun assetPath -> upAndOutEuropeanCallPayoff strike upperBarrier assetPath)

    controlVariate numTrajectories
                   discountFactor
                   (fun () -> getAssetPath S0 r deltaT sigma normal numSamples)
                   priceUpAndOutExact
                   barrierCallPayoff
                   upAndOut

let asian =
    let S0 = 100.0
    let strike = 90.0
    let r = 0.05
    let T = 1.0
    let sigma = 0.2
    let numTrajectories = 100000
    let numSamples = 12

    let result = priceAsianArithmeticMeanMC S0 strike r T sigma numTrajectories numSamples
    printfn "Asian arithmetic mean:%f stddev:%f" (fst result) (snd result)

let barrier = 
    let S0 = 100.0
    let strike = 90.0
    let upperBarrier = 160.0
    let lowerBarrier = 75.0
    let r = 0.05
    let T = 0.5
    let sigma = 0.4
    let numTrajectories = 1000
    let numSamples = 10000
    
    let result = priceDoubleBarrierCallMC S0 strike upperBarrier lowerBarrier r T sigma numTrajectories numSamples
    printfn "Monte Carlo double barrier call:%f stddev:%f" (fst result) (snd result)

    let resultMCCVEU = priceDoubleBarrierCall_MC_CV_EU S0 strike upperBarrier lowerBarrier r T sigma numTrajectories numSamples
    printfn "Control Variate with European Call price:%f stddev:%f" (fst resultMCCVEU) (snd resultMCCVEU)

    let resultMCCVUpAndOut = priceDoubleBarrierCall_MC_CV_UpAndOut S0 strike upperBarrier lowerBarrier r T sigma numTrajectories numSamples
    printfn "Control Variate with Up And Out price:%f stddev:%f" (fst resultMCCVUpAndOut) (snd resultMCCVUpAndOut)

[<EntryPoint>]
let main argv = 
    asian
    barrier
    0