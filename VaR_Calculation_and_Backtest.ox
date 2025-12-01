#include <oxstd.h>

/**
 * Variance-Covariance Approach for Value at Risk Calculation
 * Author: Generated Code
 * Description: Calculate VaR for 19 stocks using parametric method
 */

// =================== PARAMETERS (Easy to modify) ===================
decl CONFIDENCE_LEVEL;           // 99% confidence level
decl INVESTMENT_AMOUNT;          // Default investment amount per stock (currency units) - for backward compatibility
decl vINVESTMENT_AMOUNTS;        // Investment amounts for each asset (currency units vector)
decl HISTORICAL_PERIOD;          // Historical observation period (trading days)
decl HOLDING_PERIOD;             // Holding period (days)
decl DATA_FILE;                  // Price data file
decl NUM_ASSETS;                 // Number of stocks
decl USE_LOG_RETURNS;            // TRUE for log returns, FALSE for simple returns
decl VAR_METHOD;                 // "PARAMETRIC" for Variance-Covariance, "HISTORICAL" for Historical Simulation
decl g_iCallCount;               // Global counter for backtesting random seed

// Initialize parameters
initializeParameters()
{
    CONFIDENCE_LEVEL = 0.99;
    INVESTMENT_AMOUNT = 100000;  // Default amount for backward compatibility
    HISTORICAL_PERIOD = 1000;
    HOLDING_PERIOD = 10;
    DATA_FILE = "pdatau.prn";
    NUM_ASSETS = 19;
    USE_LOG_RETURNS = TRUE;      // Set to FALSE to use simple returns instead
    VAR_METHOD = "HISTORICAL";   // Set to "HISTORICAL" to use Historical Simulation
    
    // Initialize investment amounts for each asset
    // Option 1: Equal investment (default behavior)
     vINVESTMENT_AMOUNTS = ones(NUM_ASSETS, 1) * INVESTMENT_AMOUNT;
    
    // Option 2: Custom investment amounts (uncomment and modify as needed)
    // Example: Different investment amounts for each asset
    // vINVESTMENT_AMOUNTS = <0; 0; 3000000; 0; 3000000; 3000000; 3000000; 3000000; 3000000; 3000000; 3000000; 3000000; 3000000; 3000000; 0; 0; 0; 0; 3000000>;
    
    // Option 3: Market cap weighted (example)
    // vINVESTMENT_AMOUNTS = <200000; 180000; 160000; 140000; 120000;
    //                       100000; 90000; 80000; 70000; 60000;
    //                       50000; 45000; 40000; 35000; 30000;
    //                       25000; 20000; 15000; 10000>;
    
    // Initialize global counter
    g_iCallCount = 0;
}

/**
 * Load price data from file
 */
loadPriceData()
{
    decl sDataFile = "pdatau.prn";  // Use local variable to avoid initialization issues
    println("- Loading data from: ", sDataFile);
    
    decl mData = loadmat(sDataFile);
    if (mData == <>) {
        println("Error: Could not load data file: ", sDataFile);
        return <>;
    }
    
    // Skip the first row (which contains dimensions) and return price matrix
    if (rows(mData) < 2) {
        println("Error: Data file appears to be empty or invalid");
        return <>;
    }
    
    return mData[1:][]; // Skip first row with dimensions
}

/**
 * Load trading dates from timept.xls
 */
loadTradingDates()
{
    println("- Loading trading dates from: timept.xls");
    
    decl mDates = loadmat("timept.xls");
    if (mDates == <>) {
        println("Warning: Could not load date file, will use index numbers instead");
        return <>;
    }
    
    println("- Trading dates loaded: ", rows(mDates), " dates");
    return mDates;
}

/**
 * Save backtesting results with custom headers to CSV
 */
saveBacktestResultsCSV(const mResults, const mDates)
{
    decl sFileName;
    if (VAR_METHOD == "HISTORICAL") {
        sFileName = "VaR_Backtest_Results_Historical.csv";
    } else {
        sFileName = "VaR_Backtest_Results_Parametric.csv";
    }
    
    println("\nSaving detailed results to ", sFileName, "...");
    
    decl fp = fopen(sFileName, "w");
    if (fp == 0) {
        println("Error: Could not create CSV file, using default savemat");
        savemat(sFileName, mResults);
        return;
    }
    
    // Write custom headers with method information
    fprint(fp, "Test_Period,Start_Date,VaR_Predicted_", VAR_METHOD, ",Actual_Loss,Violation_Flag\n");
    
    // Write data rows
    decl i;
    for (i = 0; i < rows(mResults); i++) {
        decl iPredictionIdx = mResults[i][3];  // VaR prediction date index from results
        decl sStartDate;
        
        if (mDates != <> && iPredictionIdx < rows(mDates)) {
            // Use actual date if available
            sStartDate = sprint("%.0f", mDates[iPredictionIdx][0]);
        } else {
            // Use index if dates not available
            sStartDate = sprint("Day_%d", iPredictionIdx + 1);
        }
        
        // Write: Test_Period, Start_Date, VaR_Predicted, Actual_Loss, Violation_Flag
        fprint(fp, sprint("%d", i + 1), ",", sStartDate, ",", 
               sprint("%.6f", mResults[i][0]), ",",
               sprint("%.6f", mResults[i][1]), ",",
               sprint("%.0f", mResults[i][2]), "\n");
    }
    
    fclose(fp);
    println("Results saved with custom headers!");
    println("Columns:");
    println("- Test_Period: Sequential test number (1 to ", rows(mResults), ")");
    println("- Start_Date: Trading date for VaR prediction (start of holding period)");
    println("- VaR_Predicted_", VAR_METHOD, ": VaR predicted using ", VAR_METHOD, " method (currency units)");
    println("- Actual_Loss: Actual portfolio loss over holding period (currency units)");  
    println("- Violation_Flag: 1=violation occurred, 0=no violation");
}

/**
 * Calculate log returns from price data
 * 
 * Log Returns vs Simple Returns:
 * 
 * 1. Log Returns: r_t = ln(P_t / P_{t-1})
 *    Advantages:
 *    - Time additivity: r(t,t+k) = r(t,t+1) + r(t+1,t+2) + ... + r(t+k-1,t+k)
 *    - Better approximation to normal distribution
 *    - Symmetric: log(P_t/P_{t-1}) = -log(P_{t-1}/P_t)
 *    - More suitable for multi-period VaR scaling
 * 
 * 2. Simple Returns: r_t = (P_t - P_{t-1}) / P_{t-1}
 *    Advantages:
 *    - More intuitive interpretation (direct percentage change)
 *    - Portfolio returns are weighted averages of individual returns
 *    - Easier to understand for practitioners
 */
calculateLogReturns(const mPrices)
{
    decl iRows = rows(mPrices);
    if (iRows < 2) {
        println("Error: Need at least 2 price observations to calculate returns");
        return <>;
    }
    
    // Calculate log returns: r_t = ln(P_t / P_{t-1})
    decl mReturns = log(mPrices[1:][]) - log(mPrices[:iRows-2][]);
    
    return mReturns;
}

/**
 * Calculate simple returns from price data
 */
calculateSimpleReturns(const mPrices)
{
    decl iRows = rows(mPrices);
    if (iRows < 2) {
        println("Error: Need at least 2 price observations to calculate returns");
        return <>;
    }
    
    // Calculate simple returns: r_t = (P_t - P_{t-1}) / P_{t-1}
    decl mReturns = (mPrices[1:][]) ./ (mPrices[:iRows-2][]) - 1;
    
    return mReturns;
}

/**
 * Calculate portfolio weights based on investment amounts
 */
calculateWeights(const vPrices)
{
    // Calculate total portfolio value
    decl dTotalValue = sumc(vINVESTMENT_AMOUNTS);
    
    // Calculate weights as proportion of total investment
    decl vWeights = vINVESTMENT_AMOUNTS / dTotalValue;
    
    return vWeights;
}

/**
 * Calculate VaR using Historical Simulation method (pghs.ox style - using price differences)
 */
calculateHistoricalVaR(const mHistPrices, const vWeights)
{
    // Calculate price differences (like pghs.ox using diff0)
    decl mPriceDiffs = mHistPrices[1:][] - mHistPrices[:rows(mHistPrices)-2][];
    
    if (rows(mPriceDiffs) < HOLDING_PERIOD) {
        println("Error: Not enough historical data for simulation");
        return .NaN;
    }
    
    // Get current prices (last observation in historical data)
    decl vCurrentPrices = mHistPrices[rows(mHistPrices) - 1][];
    
    // Calculate current portfolio value using share quantities
    decl vShares = vINVESTMENT_AMOUNTS ./ vCurrentPrices;  // Number of shares for each asset
    decl dCurrentPortfolioValue = (vCurrentPrices * vShares)[0];
    
    // Monte Carlo simulation (following pghs.ox approach)
    decl iNumSimulations = 1000;  // Same as pghs.ox
    decl iNumHistoricalPeriods = rows(mPriceDiffs);
    decl vSimulatedPortfolioValues = zeros(iNumSimulations, 1);
    
    decl iSim, iDay;
    for (iSim = 0; iSim < iNumSimulations; iSim++) {
        // Start with current prices (like pghs.ox: portrl=pdt[2013+i][])
        decl vSimPrices = vCurrentPrices;
        
        // Simulate HOLDING_PERIOD days using random historical price changes
        for (iDay = 0; iDay < HOLDING_PERIOD; iDay++) {
            // Randomly select a historical period (like pghs.ox random selection)
            decl iRandomPeriod = floor(ranu(1, 1)[0] * iNumHistoricalPeriods);
            if (iRandomPeriod >= iNumHistoricalPeriods) iRandomPeriod = iNumHistoricalPeriods - 1;
            
            // Add historical price changes (like pghs.ox: +pdtdf[random_index])
            vSimPrices = vSimPrices + mPriceDiffs[iRandomPeriod][];
        }
        
        // Calculate simulated portfolio value using fixed share quantities
        vSimulatedPortfolioValues[iSim] = (vSimPrices * vShares)[0];
    }
    
    // Calculate VaR as difference between current value and 1% quantile (like pghs.ox)
    decl vSortedValues = sortc(vSimulatedPortfolioValues);
    decl iQuantileIndex = floor(0.01 * iNumSimulations);  // 1% quantile for 99% confidence
    if (iQuantileIndex >= iNumSimulations) iQuantileIndex = iNumSimulations - 1;
    
    decl dWorstCaseValue = vSortedValues[iQuantileIndex];
    decl dHistoricalVaR = dCurrentPortfolioValue - dWorstCaseValue;
    
    return dHistoricalVaR;
}

/**
 * Calculate VaR using Parametric (Variance-Covariance) method
 */
calculateParametricVaR(const mHistPrices, const vWeights)
{
    // Calculate returns
    decl mReturns;
    if (USE_LOG_RETURNS) {
        mReturns = calculateLogReturns(mHistPrices);
    } else {
        mReturns = calculateSimpleReturns(mHistPrices);
    }
    
    if (mReturns == <>) {
        return .NaN;
    }
    
    // Calculate statistics
    decl vMeanReturns = meanc(mReturns);
    decl mCovMatrix = variance(mReturns);
    
    // Ensure correct dimensions
    if (columns(vMeanReturns) > rows(vMeanReturns)) {
        vMeanReturns = vMeanReturns';
    }
    
    // Portfolio statistics
    decl dPortfolioMean = vWeights' * vMeanReturns;
    decl dPortfolioVar = vWeights' * mCovMatrix * vWeights;
    decl dPortfolioStd = sqrt(dPortfolioVar);
    
    // Scale for holding period
    decl dScaledMean = dPortfolioMean * HOLDING_PERIOD;
    decl dScaledStd = dPortfolioStd * sqrt(HOLDING_PERIOD);
    
    // Calculate VaR
    decl dAlpha = 1 - CONFIDENCE_LEVEL;
    decl dZScore = quann(dAlpha);
    decl dPortfolioValue = sumc(vINVESTMENT_AMOUNTS);
    decl dParametricVaR = -(dScaledMean + dZScore * dScaledStd) * dPortfolioValue;
    
    return dParametricVaR;
}

/**
 * Calculate VaR for a single period (used in backtesting)
 */
calculateSinglePeriodVaR(const mHistPrices)
{
    // Set random seed for consistent backtesting (increment based on period for variation)
    if (g_iCallCount == <>) {
        g_iCallCount = 0;
    }
    ranseed(967537412 + g_iCallCount);
    g_iCallCount++;
    
    // Calculate weights using latest prices
    decl vLatestPrices = mHistPrices[rows(mHistPrices) - 1][];
    decl vWeights = calculateWeights(vLatestPrices);
    
    // Choose VaR calculation method
    decl dVaR;
    if (VAR_METHOD == "HISTORICAL") {
        dVaR = calculateHistoricalVaR(mHistPrices, vWeights);
    } else {
        dVaR = calculateParametricVaR(mHistPrices, vWeights);
    }
    
    return dVaR;
}

/**
 * Calculate actual portfolio loss over a specific period
 */
calculateActualLoss(const mAllPrices, const iStart, const iEnd)
{
    // Get prices at start and end of holding period
    decl vStartPrices = mAllPrices[iStart][];
    decl vEndPrices = mAllPrices[iEnd][];
    
    // Calculate portfolio returns for each asset
    decl vReturns;
    if (USE_LOG_RETURNS) {
        vReturns = log(vEndPrices) - log(vStartPrices);
    } else {
        vReturns = (vEndPrices - vStartPrices) ./ vStartPrices;
    }
    
    // Portfolio return (weighted by investment amounts)
    decl vWeights = vINVESTMENT_AMOUNTS / sumc(vINVESTMENT_AMOUNTS);
    
    // Ensure correct dimensions: vReturns should be row vector, vWeights should be column vector
    if (rows(vReturns) > columns(vReturns)) {
        vReturns = vReturns';  // Convert to row vector if needed
    }
    if (columns(vWeights) > rows(vWeights)) {
        vWeights = vWeights';  // Convert to column vector if needed
    }
    
    decl dPortfolioReturn = (vReturns * vWeights)[0];  // Weighted portfolio return
    
    // Portfolio loss (negative of return, in currency units)
    decl dPortfolioValue = sumc(vINVESTMENT_AMOUNTS);
    decl dActualLoss = -dPortfolioReturn * dPortfolioValue;
    
    return dActualLoss;
}

/**
 * Fast version of Historical VaR calculation for backtesting (fewer simulations)
 */
calculateHistoricalVaR_Fast(const mHistPrices, const iNumSimulations)
{
    // Set random seed for consistent results
    if (g_iCallCount == <>) {
        g_iCallCount = 0;
    }
    ranseed(967537412 + g_iCallCount);
    g_iCallCount++;
    
    // Calculate price differences (like pghs.ox using diff0)
    decl mPriceDiffs = mHistPrices[1:][] - mHistPrices[:rows(mHistPrices)-2][];
    
    if (rows(mPriceDiffs) < HOLDING_PERIOD) {
        return .NaN;
    }
    
    // Get current prices (last observation in historical data)
    decl vCurrentPrices = mHistPrices[rows(mHistPrices) - 1][];
    
    // Calculate current portfolio value using share quantities
    decl vShares = vINVESTMENT_AMOUNTS ./ vCurrentPrices;
    decl dCurrentPortfolioValue = (vCurrentPrices * vShares)[0];
    
    // Fast Monte Carlo simulation
    decl iNumHistoricalPeriods = rows(mPriceDiffs);
    decl vSimulatedPortfolioValues = zeros(iNumSimulations, 1);
    
    decl iSim, iDay;
    for (iSim = 0; iSim < iNumSimulations; iSim++) {
        decl vSimPrices = vCurrentPrices;
        
        for (iDay = 0; iDay < HOLDING_PERIOD; iDay++) {
            decl iRandomPeriod = floor(ranu(1, 1)[0] * iNumHistoricalPeriods);
            if (iRandomPeriod >= iNumHistoricalPeriods) iRandomPeriod = iNumHistoricalPeriods - 1;
            vSimPrices = vSimPrices + mPriceDiffs[iRandomPeriod][];
        }
        
        vSimulatedPortfolioValues[iSim] = (vSimPrices * vShares)[0];
    }
    
    // Calculate VaR
    decl vSortedValues = sortc(vSimulatedPortfolioValues);
    decl iQuantileIndex = floor(0.01 * iNumSimulations);
    if (iQuantileIndex >= iNumSimulations) iQuantileIndex = iNumSimulations - 1;
    
    decl dWorstCaseValue = vSortedValues[iQuantileIndex];
    decl dHistoricalVaR = dCurrentPortfolioValue - dWorstCaseValue;
    
    return dHistoricalVaR;
}

/**
 * Perform VaR backtesting for all possible periods
 * Returns matrix with columns: [VaR_predicted, Actual_loss, Violation_flag, Start_date]
 */
performVaRBacktest(const mAllPrices)
{
    decl iTotalDays = rows(mAllPrices);
    decl iBacktestPeriods = iTotalDays - HISTORICAL_PERIOD - HOLDING_PERIOD;
    
    if (iBacktestPeriods <= 0) {
        println("Error: Insufficient data for backtesting");
        return <>;
    }
    
    println("Performing VaR backtesting for ", iBacktestPeriods, " periods...");
    println("Each test uses ", HISTORICAL_PERIOD, " days of history to predict ", HOLDING_PERIOD, "-day VaR");
    
    // For speed optimization, reduce number of simulations for backtesting
    decl iOriginalSimulations = 1000;
    decl iBacktestSimulations = 500;  // Reduce simulations for speed
    
    if (VAR_METHOD == "HISTORICAL") {
        println("Speed optimization: Using ", iBacktestSimulations, " simulations per period for backtesting");
        println("(Final calculation will use full ", iOriginalSimulations, " simulations)");
    }
    
    // Initialize results matrix: [VaR_predicted, Actual_loss, Violation_flag, Start_date]
    decl mResults = zeros(iBacktestPeriods, 4);
    
    decl iViolations = 0;
    decl i, iStartIdx, iEndIdx, iHoldingStart, iHoldingEnd;
    decl iProgressInterval = max(50, iBacktestPeriods / 20);  // Show progress every 5%
    
    for (i = 0; i < iBacktestPeriods; i++) {
        // Show progress at regular intervals
        if (fmod(i, iProgressInterval) == 0 || i == iBacktestPeriods - 1) {
            decl dProgress = (i + 1.0) / iBacktestPeriods * 100;
            println("Progress: ", dProgress, "% (", i + 1, "/", iBacktestPeriods, ")");
        }
        
        // Define data windows
        iStartIdx = i;
        iEndIdx = i + HISTORICAL_PERIOD - 1;
        iHoldingStart = i + HISTORICAL_PERIOD;
        iHoldingEnd = i + HISTORICAL_PERIOD + HOLDING_PERIOD - 1;
        
        // Get historical data for VaR calculation
        decl mHistPrices = mAllPrices[iStartIdx:iEndIdx][];
        
        // Calculate VaR for this period (optimized for speed during backtesting)
        decl dVaR;
        if (VAR_METHOD == "HISTORICAL") {
            dVaR = calculateHistoricalVaR_Fast(mHistPrices, iBacktestSimulations);
        } else {
            dVaR = calculateSinglePeriodVaR(mHistPrices);
        }
        
        // Calculate actual portfolio loss over holding period
        decl dActualLoss = calculateActualLoss(mAllPrices, iHoldingStart, iHoldingEnd);
        
        // Check for violation (actual loss > VaR)
        decl iViolation = (dActualLoss > dVaR) ? 1 : 0;
        if (iViolation) iViolations++;
        
        // Store results
        mResults[i][0] = dVaR;
        mResults[i][1] = dActualLoss;
        mResults[i][2] = iViolation;
        mResults[i][3] = iHoldingStart;  // Store the start of holding period (VaR prediction date)
    }
    
    println("Backtesting completed!");
    println("Total violations: ", iViolations, " out of ", iBacktestPeriods, " periods");
    println("Violation rate: ", iViolations / iBacktestPeriods * 100, " percent");
    println("Expected violation rate: ", (1 - CONFIDENCE_LEVEL) * 100, " percent");
    
    return mResults;
}



/**
 * Analyze backtesting results and provide summary statistics
 */
analyzeBacktestResults(const mResults, const mDates)
{
    decl iNumTests = rows(mResults);
    decl vVaRPredicted = mResults[][0];
    decl vActualLoss = mResults[][1];
    decl vViolations = mResults[][2];
    
    // Basic statistics
    decl iViolations = sumc(vViolations);
    decl dViolationRate = iViolations / iNumTests;
    decl dExpectedRate = 1 - CONFIDENCE_LEVEL;
    
    println("\nBacktesting Results Summary:");
    println("============================");
    println("Number of tests: ", iNumTests);
    println("Number of violations: ", iViolations);
    print("Actual violation rate: ", "%.2f", dViolationRate * 100, " percent (Expected: ", "%.2f", dExpectedRate * 100, " percent)\n");
    
    // Statistics on VaR predictions and actual losses
    println("\nVaR Predictions:");
    print("- Mean VaR: ", "%.2f", meanc(vVaRPredicted), " currency units\n");
    print("- Min VaR:  ", "%.2f", minc(vVaRPredicted), " currency units\n");
    print("- Max VaR:  ", "%.2f", maxc(vVaRPredicted), " currency units\n");
    
    println("\nActual Losses:");
    print("- Mean loss: ", "%.2f", meanc(vActualLoss), " currency units\n");
    print("- Min loss:  ", "%.2f", minc(vActualLoss), " currency units\n");
    print("- Max loss:  ", "%.2f", maxc(vActualLoss), " currency units\n");
    
    // Model performance assessment
    println("\nModel Performance Assessment:");
    if (fabs(dViolationRate - dExpectedRate) < 0.01) {
        println("✓ Model performance: GOOD - Violation rate close to expected");
    } else if (dViolationRate > dExpectedRate + 0.02) {
        println("⚠ Model performance: POOR - Too many violations (model underestimates risk)");
    } else if (dViolationRate < dExpectedRate - 0.02) {
        println("⚠ Model performance: CONSERVATIVE - Too few violations (model overestimates risk)");
    } else {
        println("△ Model performance: ACCEPTABLE - Minor deviation from expected");
    }
    
    // Save results with improved format
    saveBacktestResultsCSV(mResults, mDates);
}

// Main function - put at the end
main()
{
    // Initialize parameters first
    initializeParameters();
    
    // Set random seed for reproducible results (like pghs.ox)
    ranseed(967537412);
    
    println("========================================");
    println("Value at Risk Calculation");
    println("Method: Variance-Covariance Approach");
    println("========================================");
    
    // Display parameters
    println("Parameters:");
    println("- Confidence Level: ", CONFIDENCE_LEVEL * 100, " percent");
    println("- Total Portfolio Investment: ", sumc(vINVESTMENT_AMOUNTS), " currency units");
    println("- Investment Strategy: ", (variance(vINVESTMENT_AMOUNTS) < 0.01) ? "Equal Weight" : "Custom Weight");
    println("- Historical Period: ", HISTORICAL_PERIOD, " trading days");
    println("- Holding Period: ", HOLDING_PERIOD, " days");
    println("- Number of Assets: ", NUM_ASSETS);
    println("- Data File: ", DATA_FILE);
    println("- Return Calculation Method: ", USE_LOG_RETURNS ? "Log Returns" : "Simple Returns");
    println("- VaR Method: ", VAR_METHOD);
    println("");
    
    // Step 1: Load price data
    println("Step 1: Loading price data...");
    decl mPrices = loadPriceData();
    if (mPrices == <>) {
        println("Error: Failed to load price data!");
        return;
    }
    
    decl iRows = rows(mPrices);
    decl iCols = columns(mPrices);
    println("- Data loaded: ", iRows, " observations x ", iCols, " assets");
    
    // Step 1.5: Load trading dates
    println("\nStep 1.5: Loading trading dates...");
    decl mDates = loadTradingDates();
    
    // Check if we have enough data
    if (iRows < HISTORICAL_PERIOD + 1) {
        println("Error: Insufficient data! Need at least ", HISTORICAL_PERIOD + 1, " observations.");
        return;
    }
    
    // Step 2: Use the most recent data for calculation
    decl mRecentPrices = mPrices[iRows - HISTORICAL_PERIOD - 1 : iRows - 1][]; // Last 1001 observations
    println("- Using most recent ", HISTORICAL_PERIOD + 1, " observations for calculation");
    
    // Step 3: Calculate returns
    println("\nStep 2: Calculating returns...");
    decl mReturns;
    if (USE_LOG_RETURNS) {
        mReturns = calculateLogReturns(mRecentPrices);
        println("- Calculated ", rows(mReturns), " log return observations");
    } else {
        mReturns = calculateSimpleReturns(mRecentPrices);
        println("- Calculated ", rows(mReturns), " simple return observations");
    }
    
    if (mReturns == <>) {
        println("Error: Returns calculation failed!");
        return;
    }
    
    // Step 4: Calculate statistics
    println("\nStep 3: Calculating portfolio statistics...");
    decl vMeanReturns = meanc(mReturns);
    decl mCovMatrix = variance(mReturns);
    
    // Step 5: Calculate portfolio weights
    println("\nStep 4: Calculating portfolio weights...");
    decl vLatestPrices = mRecentPrices[rows(mRecentPrices) - 1][];
    decl vWeights = calculateWeights(vLatestPrices);
    
    println("Asset allocation details:");
    decl i;
    for (i = 0; i < NUM_ASSETS; i++) {
        decl dInvestment = vINVESTMENT_AMOUNTS[i];
        decl dHoldings = dInvestment / vLatestPrices[i];
        print("Asset ", "%2.0f", i + 1, ": Price = ", "%8.3f", vLatestPrices[i], 
              ", Investment = ", "%10.0f", dInvestment, ", Weight = ", "%7.4f", vWeights[i][0] * 100, "%",
              ", Holdings = ", "%10.3f", dHoldings, "\n");
    }
    
    // Step 6: Calculate portfolio return statistics
    println("\nStep 5: Calculating portfolio return statistics...");
    
    // Ensure correct vector dimensions for matrix multiplication
    if (columns(vMeanReturns) > rows(vMeanReturns)) {
        vMeanReturns = vMeanReturns';
    }
    if (columns(vWeights) > rows(vWeights)) {
        vWeights = vWeights';
    }
    
    decl dPortfolioMean = vWeights' * vMeanReturns;
    decl dPortfolioVar = vWeights' * mCovMatrix * vWeights;
    
    decl dPortfolioStd = sqrt(dPortfolioVar);
    
    println("Portfolio statistics (daily):");
    print("- Mean return: ", "%8.6f", dPortfolioMean, " (", "%.4f", dPortfolioMean * 100, " percent)\n");
    print("- Volatility:  ", "%8.6f", dPortfolioStd, " (", "%.4f", dPortfolioStd * 100, " percent)\n");
    
    // Step 7: Scale for holding period
    println("\nStep 6: Scaling for holding period...");
    decl dScaledMean = dPortfolioMean * HOLDING_PERIOD;
    decl dScaledStd = dPortfolioStd * sqrt(HOLDING_PERIOD);
    
    print("Portfolio statistics (", "%d", HOLDING_PERIOD, "-day holding period):\n");
    print("- Expected return: ", "%8.6f", dScaledMean, " (", "%.4f", dScaledMean * 100, " percent)\n");
    print("- Volatility:      ", "%8.6f", dScaledStd, " (", "%.4f", dScaledStd * 100, " percent)\n");
    
    // Step 8: Calculate portfolio value
    decl dPortfolioValue = sumc(vINVESTMENT_AMOUNTS);
    println("\nStep 7: Calculating Value at Risk...");
    print("- Total portfolio value: ", "%12.2f", dPortfolioValue, " currency units\n");
    
    // Step 9: Calculate VaR
    println("\nResults:");
    println("========================================");
    
    decl dVaR;
    if (VAR_METHOD == "HISTORICAL") {
        println("Using Historical Simulation method (Monte Carlo approach)...");
        dVaR = calculateHistoricalVaR(mRecentPrices, vWeights);
        println("Historical VaR calculated using 1000 Monte Carlo simulations");
        println("Each simulation samples ", HOLDING_PERIOD, " random days from ", HISTORICAL_PERIOD, " historical periods");
    } else {
        println("Using Parametric (Variance-Covariance) method...");
        decl dAlpha = 1 - CONFIDENCE_LEVEL;
        decl dZScore = quann(dAlpha);  // Critical value from normal distribution
        dVaR = -(dScaledMean + dZScore * dScaledStd) * dPortfolioValue;
        print("Critical z-score (", "%.1f", CONFIDENCE_LEVEL * 100, " percent confidence): ", "%8.4f", dZScore, "\n");
    }
    
    decl dVaRPercentage = dVaR / dPortfolioValue * 100;
    
    println("Value at Risk (VaR):");
    print("- Absolute amount: ", "%12.2f", dVaR, " currency units\n");
    print("- Percentage of portfolio: ", "%8.4f", dVaRPercentage, " percent\n");
    print("- Interpretation: With ", "%.1f", CONFIDENCE_LEVEL * 100, 
          " percent confidence, the portfolio will not lose more than ", "%.2f", dVaR, 
          " currency units over ", "%d", HOLDING_PERIOD, " days\n");
    
    // Step 9.5: VaR Method Comparison
    println("\n========================================");
    println("VaR Method Comparison");
    println("========================================");
    
    decl dParametricVaR = calculateParametricVaR(mRecentPrices, vWeights);
    decl dHistoricalVaR = calculateHistoricalVaR(mRecentPrices, vWeights);
    decl dDifference = dHistoricalVaR - dParametricVaR;
    decl dDifferencePercent = dDifference / dParametricVaR * 100;
    
    print("Parametric VaR (Variance-Covariance): ", "%12.2f", dParametricVaR, " currency units\n");
    print("Historical VaR (Historical Simulation): ", "%12.2f", dHistoricalVaR, " currency units\n");
    print("Difference (Historical - Parametric): ", "%12.2f", dDifference, " currency units\n");
    print("Relative Difference: ", "%8.2f", dDifferencePercent, " percent\n");
    
    if (fabs(dDifferencePercent) < 5) {
        println("→ Methods show similar results (difference < 5 percent)");
    } else if (dHistoricalVaR > dParametricVaR) {
        println("→ Historical method suggests higher risk than parametric method");
    } else {
        println("→ Historical method suggests lower risk than parametric method");
    }
    
    // Additional risk metrics
    println("\nAdditional Risk Metrics:");
    print("- Daily VaR (1-day): ", "%12.2f", dVaR / sqrt(HOLDING_PERIOD), " currency units\n");
    print("- VaR as percent of total investment: ", "%8.4f", dVaR / dPortfolioValue * 100, " percent\n");
    
    println("\n========================================");
    println("Calculation completed successfully!");
    println("\nNote: To modify investment weights, edit the vINVESTMENT_AMOUNTS");
    println("      vector in the initializeParameters() function.");
    println("      Current setup: ", (variance(vINVESTMENT_AMOUNTS) < 0.01) ? "Equal weights" : "Custom weights");
    
    // Step 10: VaR Backtesting
    println("\n\n========================================");
    println("VaR Backtesting Analysis");
    println("========================================");
    
    decl mBacktestResults = performVaRBacktest(mPrices);
    if (mBacktestResults != <>) {
        analyzeBacktestResults(mBacktestResults, mDates);
    }
} 