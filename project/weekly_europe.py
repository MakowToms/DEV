import json
from icecream import ic
import pandas as pd
from natsort import natsorted
import matplotlib.pyplot as plt
import datetime

data = json.load(open('project/data/weekly_surveillance_stats.json', 'r'))
data = data['stats']

keys = ['submissions_per_variant', 'submissions_per_lineage', 'submissions_per_clade', 'submissions']

# number of submissions in weeks for top7 countries
result = {}

for week in data.keys():
    week_data = data[week]
    week_res = {}
    for country in week_data.keys():
        country_data = week_data[country]
        week_res[country] = country_data[keys[3]]
    week_name = datetime.datetime.strptime(week + '-1', "%W-%Y-%w")
    result[week_name] = week_res

df = pd.DataFrame(result)
df = df[natsorted(df.columns)]
df = df.fillna(0)
df = df.transpose()

sums = df.sum()
sums = pd.DataFrame(sums)
sums = sums.sort_values(0, ascending=False)
countries = sums[:20].index
index = df.iloc[4:, :].index

df.iloc[4:, :].loc[:, countries].plot()
plt.show()


# percentage of each clade in top 7 countries

percentage_res = {}
countries = list(countries)
countries.append('Europe')
countries = pd.Index(countries)
for country in countries:
    percentage_res[country] = {}

for week in data.keys():
    week_data = data[week]
    week_name = datetime.datetime.strptime(week + '-1', "%W-%Y-%w")
    percentage_res['Europe'][week_name] = {'G': 0, 'GH': 0, 'GR': 0, 'GRY': 0, 'GV': 0,
                                           'L': 0, 'O': 0, 'S': 0, 'V': 0, None: 0}
    Europe_sum = 0
    for country in week_data.keys():
        if country in countries:
            country_data = week_data[country]
            week_res = {'G': 0, 'GH': 0, 'GR': 0, 'GRY': 0, 'GV': 0,
                        'L': 0, 'O': 0, 'S': 0, 'V': 0, None: 0}
            for clade in country_data[keys[2]]:
                if not week_res.__contains__(clade['value']):
                    print(clade['value'])
                week_res[clade['value']] = clade['count'] / country_data[keys[3]]
                percentage_res['Europe'][week_name][clade['value']] += clade['count']
                Europe_sum += clade['count']
            percentage_res[country][week_name] = week_res

    for clade in percentage_res['Europe'][week_name].keys():
        percentage_res['Europe'][week_name][clade] /= Europe_sum

fig, axes = plt.subplots(countries.shape[0], sharex=True, figsize=(10, 20))
plt.subplots_adjust(hspace=0.8)  # bottom=0.9
for i, country in enumerate(countries):
    for week_name in index:
        if not percentage_res[country].__contains__(week_name):
            percentage_res[country][week_name] = {}
    df = pd.DataFrame(percentage_res[country])
    df = df[natsorted(df.columns)]
    df = df.loc[natsorted(df.index), :]
    df = df.fillna(0)
    df = df.transpose()
    df = df.loc[index, :]
    df.index = [i.date() for i in index]

    ax = df.plot(kind='bar', stacked=True,
                               title=f'Country {country}',
                               ax=axes[i], legend=i == 6)
plt.show()
