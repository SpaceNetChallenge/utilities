import pandas as pd
from pandas.stats.api import ols

filename = '/Users/dlindenbaum/cosmiQGit/visualizer-2.0/data/resultsEvaluateScene_Vegas650_iter26000.csv'

testResults = pd.read_csv(filename, header=66)
testSceneSummary = pd.read_csv('/Users/dlindenbaum/cosmiQGit/visualizer-2.0/spacenetV2_Test_SceneSummary.csv')


testResults.set_index(['ImageId'], inplace=True)
testSceneSummary.set_index(['ImageId'], inplace=True)
plotResults = pd.concat([testSceneSummary, testResults], axis=1)


# Filter by AOI

vegasIndex = plotResults.index.str.contains('Vegas')
shanghaiIndex = plotResults.index.str.contains('Shanghai')
khartoumIndex = plotResults.index.str.contains('Khartoum')
parisIndex = plotResults.index.str.contains('Paris')



plotResults[vegasIndex].plot.scatter('F1Score', 'Count')
testOLS = ols(y=plotResults[vegasIndex]['F1Score'], x=plotResults[vegasIndex]['Count'])
