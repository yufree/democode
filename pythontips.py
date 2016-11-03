# basic

# This is a comment!
# Storing a string to the variable s
s = "Hello world, world"
# printing the type. It shows 'str' which means it is a string type
print type(s)
# printing the number of characters in the string.
print len(s)
# Python gave me an error below! I guess we can't change individual parts of a string this way.
s[0] = 'h'
# from here continue writing a comment for each line explaining what the following line does.
s2 = s.replace("world", "python")
s3 = s2.replace("Hello","monty")
print s
print s2
print s3
print s[6:11]
print s[6:]
print s[-2:]
s4 = s + ' ' + s3
print s4
print s4.find('world')
print 'A string with value {0} and {1}'.format(10,20.3)
help(str)

# list

values = ['1',2,3.0,False]
print len(values)
print type(values)
print values
print values[1]
print values[:3]
print values[2:]
l = []
l.append(8)
l.append(10)
l.append(10)
l.append(12)
print l
l.remove(10)
print l
l.remove(l[0])
print l
l = range(0,10,2)
print l
l = range(-5,5)
print l
line = "This is a    \t list -      \t of strings"
print len(line.split('\t'))
print line.split('\t')
print ['add another field' for i in range(10)]
l = range(-5,5)
print l
print l.sort(reverse=True)
print l

# Tuples

t = (10,40.0,"A")
print type(t), len(t)
t[1] = 'B'

# dictinaries

data = {}
print data
print len(data), type(data)
data['k4'] = 100
data = {'number': 10, 1:'string'}
data['c'] = [1,2,3,4]
print data
print data[1]
print data['c'][3]
print data['number']

# class

class Car():
    def __init__(self, model='Ford'):
        self.model = model
        self.running = False
    def start(self):
        if self.running != True:
            print 'The car started!'
            self.running = True
        else:
            print 'The car is already running!'
    def stop(self):
        if self.running == True:
            print 'The car stopped!'
            self.running = False
        else:
            print 'The car was not running!'

ford = Car()
nissan = Car(model = 'Nissan')
ford.running
ford.start()
ford.running
nissan.running
nissan.stop()

# import data from web
# mkdir nytimes_data
# cd nytimes_data
# curl -s 'http://api.nytimes.com/svc/search/v2/articlesearch.json?q=malaysia&facet_field=source&begin_date=20120101&end_date=20121231&facet_filter=true&api-key=24120b778a5da6ae998b8171dfde375e:17:70597616' > malaysia_articles.json
# 
# cat malaysia_articles.json | python -mjson.tool >> malaysia_pretty_printed.json
# 
# cat malaysia_pretty_printed.json | grep 'pub_date'
# 
# cat malaysia_pretty_printed.json | grep 'pub_date' | sort | uniq -c

#! /usr/bin/env python
import json

nytimes = json.load(open('malaysia_pretty_printed.json', 'r'))
print ','.join(['contributor', 'pub_date', 'section', 'subsection', 'word_count', 'url'])
for article in nytimes['response']['docs']:
    print '"' + '","'.join(map(str,[article['byline'],
        article['pub_date'],
        article['section_name'],
        article['subsection_name'],
        article['word_count'],
        article['web_url']])) + '"'

# chmod +x ./nytimes_parser.py
# ./nytimes_parser.py > nytimes.csv
# less nytimes.csv

import json
import urllib2

class NYTimesScraper():
    def __init__(self, apikey):
        # Creates a new NYTimesScraper Object using the apikey that was included.
        self.key = apikey
        self.url = 'http://api.nytimes.com/svc/search/v2/articlesearch.json?'
    def _build_params(self, params):
        if not params:
            raise Exception('no search parameters!')
        else:
            return '&'.join([k + '=' + v for k,v in params.iteritems()])
    def search(self, params={}):
        url = self.url + self._build_params(params)
        url = url + '&api-key=%s' % self.key
        req = urllib2.Request(url)
        data = urllib2.urlopen(req).read()
        return json.loads(data)

nytimes = NYTimesScraper(apikey='24120b778a5da6ae998b8171dfde375e:17:70597616')
articles = nytimes.search({'q':'malaysia', 'begin_date': '20140101'})

filename = 'nytimesdata.csv'
data = open(filename, 'w')
for article in articles['response']['docs']:       
    data.write('"' + '","'.join(map(str,[article['byline']['person'][0]['lastname'],
    article['pub_date'],
    article['section_name'],
    article['subsection_name'],
    article['word_count'],
    article['web_url']])) + '"\n')

data.close()
